
#include <Python.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <cvodes/cvodes.h>
#include <nvector/nvector_serial.h>
#define MYOKIT_SUNDIALS_VERSION 30100
#if MYOKIT_SUNDIALS_VERSION >= 30000
    #include <sunmatrix/sunmatrix_dense.h>
    #include <sunlinsol/sunlinsol_dense.h>
    #include <cvodes/cvodes_direct.h>
#else
    #include <cvodes/cvodes_dense.h>
#endif
#include <sundials/sundials_types.h>

#include "pacing.h"










/*

#TODO CHECK IF STILL TRUE/RELEVANT
#
# About the bindings:
#
# Time is bound to "t", the time variable used in the function. This is
#  required for periodic logging and point-list logging, when "tlog" increases
#  while engine_time stays fixed at the current solver time.
# Pace is bound to engine_pace, since the solver always visits the points where
#  its value changes the same problem as with logging "engine_time" doesn't
#  occur.
# Realtime is only useful without variable logging, so again binding to global
#  is ok.
# Evaluations may increase during interpolation, but this evaluation will be
#  taken into account in the global variable "engine_evaluations", so this is
#  fine.
#

*/





/*
 * Model error flags
 */
typedef int Model_Flag;
#define Model_OK                             0
#define Model_OUT_OF_MEMORY                 -1
// General
#define Model_INVALID_MODEL                 -10
// ESys_Populate
// ESys_AdvanceTime
// ESys_ScheduleEvent

/*
 * Sets a python exception based on a model flag.
 *
 * Arguments
 *  flag : The model flag to base the message on.
 */
void
Model_SetPyErr(Model_Flag flag)
{
    switch(flag) {
    case Model_OK:
        break;
    case Model_OUT_OF_MEMORY:
        PyErr_SetString(PyExc_Exception, "Model error: Memory allocation failed.");
        break;
    // General
    case Model_INVALID_MODEL:
        PyErr_SetString(PyExc_Exception, "Model error: Invalid model pointer provided.");
        break;
    // Populate
    // Unknown
    default:
    {
        int i = (int)flag;
        char buffer[1024];
        sprintf(buffer, "Model error: Unlisted error %d", i);
        PyErr_SetString(PyExc_Exception, buffer);
        break;
    }};
}

/*
 * Model object.
 *
 * A model struct contains all of a model's variables, including states, bound
 * variables, derivatives, and sensitivities.
 *
 * To calculate the model's RHS, the requested state and bound variables are
 * copied into the model, the RHS calculation method is called, and then the
 * calculated derivatives can be read off.
 *
 *
 */
struct Model_Mem {
    /* Number of states */
    int n_states;

    /* Whether or not this model is an ODE */
    int is_ode;

    /* Number of parameters and initial states to calculate sensitivities w.r.t. */
    int n_sens;

    /* Number of initial values to calculate sensitivities w.r.t. */
    int n_initials;

    /* Number of parameters to calculate sensitivities w.r.t. */
    int n_parameters;

    /* Number of literal constants */
    int n_literals;

    /* Bound variables */
    realtype V_time;
    realtype V_pace;

    /* Literal-derived constants */
    realtype C_RTF;
    realtype C_gK;
    realtype C_ik_IK_E;
    realtype C_ENa;
    realtype C_xp1;
    realtype C_ik1_E;
    realtype C_gK1;

    /* Parameter-derived constants */

    /* Intermediary variables */
    realtype V_i_ion;
    realtype V_i_stim;
    realtype V_ik_x_alpha;
    realtype V_ik_x_beta;
    realtype V_xi;
    realtype V_IK;
    realtype V_a;
    realtype V_ina_m_alpha;
    realtype V_ina_m_beta;
    realtype V_ina_h_alpha;
    realtype V_ina_h_beta;
    realtype V_ina_j_alpha;
    realtype V_ina_j_beta;
    realtype V_INa;
    realtype V_Kp;
    realtype V_IKp;
    realtype V_ica_E;
    realtype V_ica_d_alpha;
    realtype V_ica_d_beta;
    realtype V_ica_f_alpha;
    realtype V_ica_f_beta;
    realtype V_ICa;
    realtype V_g;
    realtype V_ik1_g_alpha;
    realtype V_ik1_g_beta;
    realtype V_IK1;
    realtype V_Ib;

    /* Sensitivity output variables */
    realtype S0_V_INa;
    realtype S1_V_INa;
    realtype S2_V_INa;

    /* State variables and derivatives */
    realtype* states;
    realtype* derivatives;

    /* Sensitivity of state variables and derivatives
     * dy/dx, where x is a parameter or initial value */
    realtype* s_states;
    realtype* s_derivatives;

    /* State-indices of states used as initial value in sensitivity analysis */
    int* i_initials;

    /* Sensitivity parameter variables */
    realtype* parameters;

    /* Literal constants */
    realtype* literals;
};
typedef struct Model_Mem* Model;

/*
 * State-variable aliases
 */
#define Y_V states[0]
#define Y_m states[1]
#define Y_h states[2]
#define Y_j states[3]
#define Y_d states[4]
#define Y_f states[5]
#define Y_x states[6]
#define Y_Ca_i states[7]
#define D_V derivatives[0]
#define D_m derivatives[1]
#define D_h derivatives[2]
#define D_j derivatives[3]
#define D_d derivatives[4]
#define D_f derivatives[5]
#define D_x derivatives[6]
#define D_Ca_i derivatives[7]

/*
 * State-variable sensitivity aliases
 */
#define S0_Y_V s_states[0]
#define S0_Y_m s_states[1]
#define S0_Y_h s_states[2]
#define S0_Y_j s_states[3]
#define S0_Y_d s_states[4]
#define S0_Y_f s_states[5]
#define S0_Y_x s_states[6]
#define S0_Y_Ca_i s_states[7]
#define S0_D_V s_derivatives[0]
#define S0_D_m s_derivatives[1]
#define S0_D_h s_derivatives[2]
#define S0_D_j s_derivatives[3]
#define S0_D_d s_derivatives[4]
#define S0_D_f s_derivatives[5]
#define S0_D_x s_derivatives[6]
#define S0_D_Ca_i s_derivatives[7]

#define S1_Y_V s_states[8]
#define S1_Y_m s_states[9]
#define S1_Y_h s_states[10]
#define S1_Y_j s_states[11]
#define S1_Y_d s_states[12]
#define S1_Y_f s_states[13]
#define S1_Y_x s_states[14]
#define S1_Y_Ca_i s_states[15]
#define S1_D_V s_derivatives[8]
#define S1_D_m s_derivatives[9]
#define S1_D_h s_derivatives[10]
#define S1_D_j s_derivatives[11]
#define S1_D_d s_derivatives[12]
#define S1_D_f s_derivatives[13]
#define S1_D_x s_derivatives[14]
#define S1_D_Ca_i s_derivatives[15]

#define S2_Y_V s_states[16]
#define S2_Y_m s_states[17]
#define S2_Y_h s_states[18]
#define S2_Y_j s_states[19]
#define S2_Y_d s_states[20]
#define S2_Y_f s_states[21]
#define S2_Y_x s_states[22]
#define S2_Y_Ca_i s_states[23]
#define S2_D_V s_derivatives[16]
#define S2_D_m s_derivatives[17]
#define S2_D_h s_derivatives[18]
#define S2_D_j s_derivatives[19]
#define S2_D_d s_derivatives[20]
#define S2_D_f s_derivatives[21]
#define S2_D_x s_derivatives[22]
#define S2_D_Ca_i s_derivatives[23]


/*
 * Parameter aliases
 */
#define P_gNa parameters[0]
#define P_gCa parameters[1]

/*
 * Literal aliases
 */
#define C_Eb literals[0]
#define C_gb literals[1]
#define C_Ca_o literals[2]
#define C_F literals[3]
#define C_K_i literals[4]
#define C_K_o literals[5]
#define C_Na_i literals[6]
#define C_Na_o literals[7]
#define C_R literals[8]
#define C_T literals[9]
#define C_PNa_K literals[10]
#define C_p1 literals[11]
#define C_gKp literals[12]
#define C_C literals[13]
#define C_i_diff literals[14]
#define C_stim_amplitude literals[15]

/*
 * Creates and returns a model struct.
 *
 * Arguments
 *  flag : The address of a model flag or NULL.
 *
 * Returns a Model pointer.
 */
Model Model_Create(Model_Flag* flag)
{
    /* Allocate */
    Model model = (Model)malloc(sizeof(struct Model_Mem));
    if (model == 0) {
        if (flag != 0) { *flag = Model_OUT_OF_MEMORY; }
        return 0;
    }

    /* Number of states */
    model->n_states = 8;
    printf("\nCreating\n");
    printf("Set n_states: %d\n", model->n_states);

    /* Number of literal constants */
    model->n_literals = 16;

    /* Whether or not this model is an ODE */
    model->is_ode = 1;

    /* Number of parameters and initial states to calculate sensitivities w.r.t. */
    model->n_sens = 3;

    /* Number of parameters to calculate sensitivities w.r.t. */
    model->n_parameters = 2;

    /* Number of initial values to calculate sensitivities w.r.t. */
    model->n_initials = 1;

    /* Set default values of bound variables */
    model->V_time = 0;
    model->V_pace = 0;

    /* State variables and derivatives */
    model->states = NULL;
    model->derivatives = NULL;
    if (model->n_states) {
        model->states = (realtype*)malloc(model->n_states * sizeof(realtype));
        model->derivatives = (realtype*)malloc(model->n_states * sizeof(realtype));
    }

    /* Sensitivity of state variables and derivatives */
    model->s_states = NULL;
    model->s_derivatives = NULL;
    if (model->n_sens && model->n_states) {
        model->s_states = (realtype*)malloc(model->n_states * model->n_sens * sizeof(realtype));
        model->s_derivatives = (realtype*)malloc(model->n_states * model->n_sens * sizeof(realtype));
    }

    /* Sensitivity parameters */
    model->parameters = NULL;
    if (model->n_parameters) {
        model->parameters = (realtype*)malloc(model->n_parameters * sizeof(realtype));
    }

    /* State-indices of states used as initial value in sensitivity analysis */
    model->i_initials = NULL;
    if (model->n_initials) {
        model->i_initials = (int*)malloc(model->n_initials * sizeof(int));
    }
    model->i_initials[0] = 7;


    /* Literal constants */
    model->literals = NULL;
    if (model->n_literals) {
        model->literals = (realtype*)malloc(model->n_literals * sizeof(realtype));
    }

    /* Set flag to indicate success */
    if (flag != 0) { *flag = Model_OK; }



    /* Return newly created model */
    return model;
}

/*
 *
 *
 * Arguments
 *  model : The model to destroy.
 *
 * Returns a model flag.
 */
Model_Flag
Model_Destroy(Model model)
{
    if (model == NULL) return Model_INVALID_MODEL;
    if (model->states != NULL) {
        free(model->states);
        model->states = NULL;
    }
    if (model->derivatives != NULL) {
        free(model->derivatives);
        model->derivatives = NULL;
    }
    if (model->parameters != NULL) {
        free(model->parameters);
        model->parameters = NULL;
    }
    if (model->i_initials != NULL) {
        free(model->i_initials);
        model->i_initials = NULL;
    }
    free(model);
    return Model_OK;
}

/*
 * Sets this model's bound variables to the given values.
 *
 * Arguments
 *  model : The model to update
 *  time
 *  pace
 *  realtime
 *  evaluations
 *
 * Returns a model flag.
 */
Model_Flag
Model_UpdateBoundVariables(
    Model model,
    realtype time, realtype pace, realtype realtime, realtype evaluations)
{
    if (model == NULL) return Model_INVALID_MODEL;

    model->V_time = time;
    model->V_pace = pace;

    return Model_OK;
}

/*
 * (Re)calculates the values of all constants that are derived from other
 * constants.
 *
 * Arguments
 *  model : The model to update.
 *
 * Returns a model flag.
 */
Model_Flag
Model_UpdateLiteralDerivedVariables(Model model)
{
    if (model == NULL) return Model_INVALID_MODEL;
    model->C_RTF = model->C_R * model->C_T / model->C_F;
    model->C_gK = 0.282 * sqrt(model->C_K_o / 5.4);
    model->C_ik_IK_E = model->C_RTF * log((model->C_K_o + model->C_PNa_K * model->C_Na_o) / (model->C_K_i + model->C_PNa_K * model->C_Na_i));
    model->C_ENa = model->C_RTF * log(model->C_Na_o / model->C_Na_i);
    model->C_xp1 = 2.0 + model->C_p1;
    model->C_ik1_E = model->C_RTF * log(model->C_K_o / model->C_K_i);
    model->C_gK1 = 0.6047 * sqrt(model->C_K_o / 5.4);

    return Model_OK;
}

/*
 * (Re)calculates the values of all constants that depend on parameters.
 *
 * Arguments
 *  model : The model to update.
 *
 * Returns a model flag.
 */
Model_Flag
Model_UpdateParameterDerivedVariables(Model model)
{
    if (model == NULL) return Model_INVALID_MODEL;

    return Model_OK;
}

/*
 * (Re)calculates the values of all intermediary variables and state
 * derivatives.
 *
 * Arguments
 *  model : The model to update
 *
 * Returns a model flag.
 */
Model_Flag
Model_UpdateDerivatives(Model model)
{
    if (model == NULL) return Model_INVALID_MODEL;

    /* ib */
    model->V_Ib = model->C_gb * (model->Y_V - model->C_Eb);
    
    /* ik */
    model->V_xi = ((model->Y_V < (-100.0)) ? 1.0 : ((model->Y_V == (-77.0)) ? 2.837 * 0.04 / exp(0.04 * (model->Y_V + 35.0)) : 2.837 * (exp(0.04 * (model->Y_V + 77.0)) - 1.0) / ((model->Y_V + 77.0) * exp(0.04 * (model->Y_V + 35.0)))));
    model->V_ik_x_alpha = 0.0005 * exp(0.083 * (model->Y_V + 50.0)) / (1.0 + exp(0.057 * (model->Y_V + 50.0)));
    model->V_ik_x_beta = 0.0013 * exp((-0.06) * (model->Y_V + 20.0)) / (1.0 + exp((-0.04) * (model->Y_V + 20.0)));
    model->D_x = model->V_ik_x_alpha * (1.0 - model->Y_x) - model->V_ik_x_beta * model->Y_x;
    model->V_IK = model->C_gK * model->V_xi * model->Y_x * (model->Y_V - model->C_ik_IK_E);
    
    /* ina */
    model->V_a = 1.0 - 1.0 / (1.0 + exp((-(model->Y_V + 40.0)) / 0.24));
    model->V_ina_m_alpha = 0.32 * (model->Y_V + 47.13) / (1.0 - exp((-0.1) * (model->Y_V + 47.13)));
    model->V_ina_m_beta = 0.08 * exp((-model->Y_V) / 11.0);
    model->D_m = model->V_ina_m_alpha * (1.0 - model->Y_m) - model->V_ina_m_beta * model->Y_m;
    model->V_INa = model->P_gNa * pow(model->Y_m, 3.0) * model->Y_h * model->Y_j * (model->Y_V - model->C_ENa);
    model->V_ina_j_alpha = model->V_a * ((-127140.0) * exp(0.2444 * model->Y_V) - 3.474e-05 * exp((-0.04391) * model->Y_V)) * (model->Y_V + 37.78) / (1.0 + exp(0.311 * (model->Y_V + 79.23)));
    model->V_ina_j_beta = model->V_a * (0.1212 * exp((-0.01052) * model->Y_V) / (1.0 + exp((-0.1378) * (model->Y_V + 40.14)))) + (1.0 - model->V_a) * (0.3 * exp((-2.535e-07) * model->Y_V) / (1.0 + exp((-0.1) * (model->Y_V + 32.0))));
    model->D_j = model->V_ina_j_alpha * (1.0 - model->Y_j) - model->V_ina_j_beta * model->Y_j;
    model->V_ina_h_alpha = model->V_a * 0.135 * exp((80.0 + model->Y_V) / (-6.8));
    model->V_ina_h_beta = model->V_a * (model->C_xp1 * exp(0.079 * model->Y_V) + 310000.0 * exp(0.35 * model->Y_V)) + (1.0 - model->V_a) / (0.13 * (1.0 + exp((model->Y_V + 10.66) / (-11.1))));
    model->D_h = model->V_ina_h_alpha * (1.0 - model->Y_h) - model->V_ina_h_beta * model->Y_h;
    
    /* ica */
    model->V_ica_E = 7.7 - 13.0287 * log(model->Y_Ca_i / model->C_Ca_o);
    model->V_ica_d_alpha = 0.095 * exp((-0.01) * (model->Y_V - 5.0)) / (1.0 + exp((-0.072) * (model->Y_V - 5.0)));
    model->V_ica_d_beta = 0.07 * exp((-0.017) * (model->Y_V + 44.0)) / (1.0 + exp(0.05 * (model->Y_V + 44.0)));
    model->D_d = model->V_ica_d_alpha * (1.0 - model->Y_d) - model->V_ica_d_beta * model->Y_d;
    model->V_ica_f_alpha = 0.012 * exp((-0.008) * (model->Y_V + 28.0)) / (1.0 + exp(0.15 * (model->Y_V + 28.0)));
    model->V_ica_f_beta = 0.0065 * exp((-0.02) * (model->Y_V + 30.0)) / (1.0 + exp((-0.2) * (model->Y_V + 30.0)));
    model->D_f = model->V_ica_f_alpha * (1.0 - model->Y_f) - model->V_ica_f_beta * model->Y_f;
    model->V_ICa = model->P_gCa * model->Y_d * model->Y_f * (model->Y_V - model->V_ica_E);
    model->D_Ca_i = (-0.0001) * model->V_ICa + 0.07 * (0.0001 - model->Y_Ca_i);
    
    /* ik1 */
    model->V_ik1_g_alpha = 1.02 / (1.0 + exp(0.2385 * (model->Y_V - model->C_ik1_E - 59.215)));
    model->V_ik1_g_beta = (0.49124 * exp(0.08032 * (model->Y_V - model->C_ik1_E + 5.476)) + 1.0 * exp(0.06175 * (model->Y_V - model->C_ik1_E - 594.31))) / (1.0 + exp((-0.5143) * (model->Y_V - model->C_ik1_E + 4.753)));
    model->V_g = model->V_ik1_g_alpha / (model->V_ik1_g_alpha + model->V_ik1_g_beta);
    model->V_IK1 = model->C_gK1 * model->V_g * (model->Y_V - model->C_ik1_E);
    
    /* ikp */
    model->V_Kp = 1.0 / (1.0 + exp((7.488 - model->Y_V) / 5.98));
    model->V_IKp = model->C_gKp * model->V_Kp * (model->Y_V - model->C_ik1_E);
    
    /* membrane */
    model->V_i_ion = model->V_INa + model->V_IK + model->V_Ib + model->V_IKp + model->V_IK1 + model->V_ICa;
    model->V_i_stim = model->V_pace * model->C_stim_amplitude;
    model->D_V = (-(1.0 / model->C_C)) * (model->V_i_ion + model->C_i_diff + model->V_i_stim);
    

    return Model_OK;
}












/*
 * Define type for "user data" that will hold parameter values if doing
 * sensitivity analysis.
 */
typedef struct {
    realtype *p;
} *UserData;

/*
 * Check sundials flags, set python error.
 *  flagvalue : The value to check
 *  funcname : The name of the function that returned the flag
 *  opt : Mode selector
 *         0 : Error if the flag is null
 *         1 : Error if the flag is < 0
 *         2 : Errir
 */
static int
check_cvode_flag(void *flagvalue, char *funcname, int opt)
{
    if (opt == 0 && flagvalue == NULL) {
        /* Check if sundials function returned null pointer */
        char str[200];
        sprintf(str, "%s() failed - returned NULL pointer", funcname);
        PyErr_SetString(PyExc_Exception, str);
        return 1;
    } else if (opt == 1) {
        /* Check if flag < 0 */
        int flag = *((int*)flagvalue);
        if (flag < 0) {
            if (strcmp(funcname, "CVode") == 0) {
                switch (flag) {
                case -1:
                    PyErr_SetString(PyExc_Exception, "Function CVode() failed with flag -1 CV_TOO_MUCH_WORK: The solver took mxstep internal steps but could not reach tout.");
                    break;
                case -2:
                    PyErr_SetString(PyExc_Exception, "Function CVode() failed with flag -2 CV_TOO_MUCH_ACC: The solver could not satisfy the accuracy demanded by the user for some internal step.");
                    break;
                case -3:
                    PyErr_SetString(PyExc_Exception, "Function CVode() failed with flag -3 CV_ERR_FAILURE: Error test failures occurred too many times during one internal time step or minimum step size was reached.");
                    break;
                case -4:
                    PyErr_SetString(PyExc_ArithmeticError, "Function CVode() failed with flag -4 CV_CONV_FAILURE: Convergence test failures occurred too many times during one internal time step or minimum step size was reached.");
                    break;
                case -5:
                    PyErr_SetString(PyExc_Exception, "Function CVode() failed with flag -5 CV_LINIT_FAIL: The linear solver's initialization function failed.");
                    break;
                case -6:
                    PyErr_SetString(PyExc_Exception, "Function CVode() failed with flag -6 CV_LSETUP_FAIL: The linear solver's setup function failed in an unrecoverable manner.");
                    break;
                case -7:
                    PyErr_SetString(PyExc_Exception, "Function CVode() failed with flag -7 CV_LSOLVE_FAIL: The linear solver's solve function failed in an unrecoverable manner.");
                    break;
                case -8:
                    PyErr_SetString(PyExc_Exception, "Function CVode() failed with flag -8 CV_RHSFUNC_FAIL: The right-hand side function failed in an unrecoverable manner.");
                    break;
                case -9:
                    PyErr_SetString(PyExc_Exception, "Function CVode() failed with flag -9 CV_FIRST_RHSFUNC_ERR: The right-hand side function failed at the first call.");
                    break;
                case -10:
                    PyErr_SetString(PyExc_Exception, "Function CVode() failed with flag -10 CV_REPTD_RHSFUNC_ERR: The right-hand side function had repeated recoverable errors.");
                    break;
                case -11:
                    PyErr_SetString(PyExc_Exception, "Function CVode() failed with flag -11 CV_UNREC_RHSFUNC_ERR: The right-hand side function had a recoverable error, but no recovery is possible.");
                    break;
                case -12:
                    PyErr_SetString(PyExc_Exception, "Function CVode() failed with flag -12 CV_RTFUNC_FAIL: The rootfinding function failed in an unrecoverable manner.");
                    break;
                case -20:
                    PyErr_SetString(PyExc_Exception, "Function CVode() failed with flag -20 CV_MEM_FAIL: A memory allocation failed.");
                    break;
                case -21:
                    PyErr_SetString(PyExc_Exception, "Function CVode() failed with flag -21 CV_MEM_NULL: The cvode mem argument was NULL.");
                    break;
                case -22:
                    PyErr_SetString(PyExc_Exception, "Function CVode() failed with flag -22 CV_ILL_INPUT: One of the function inputs is illegal.");
                    break;
                case -23:
                    PyErr_SetString(PyExc_Exception, "Function CVode() failed with flag -23 CV_NO_MALLOC: The cvode memory block was not allocated by a call to CVodeMalloc.");
                    break;
                case -24:
                    PyErr_SetString(PyExc_Exception, "Function CVode() failed with flag -24 CV_BAD_K: The derivative order k is larger than the order used.");
                    break;
                case -25:
                    PyErr_SetString(PyExc_Exception, "Function CVode() failed with flag -25 CV_BAD_T: The time t is outside the last step taken.");
                    break;
                case -26:
                    PyErr_SetString(PyExc_Exception, "Function CVode() failed with flag -26 CV_BAD_DKY: The output derivative vector is NULL.");
                    break;
                case -27:
                    PyErr_SetString(PyExc_Exception, "Function CVode() failed with flag -27 CV_TOO_CLOSE: The output and initial times are too close to each other.");
                    break;
                default: {
                     /* Note: Brackets are required here, default: should be followed by
                        a _statement_ and char str[200]; is technically not a statement... */
                    char str[200];
                    sprintf(str, "Function CVode() failed with unknown flag = %d", flag);
                    PyErr_SetString(PyExc_Exception, str);
                }}
            } else {
                char str[200];
                sprintf(str, "%s() failed with flag = %d", funcname, flag);
                PyErr_SetString(PyExc_Exception, str);
            }
            return 1;
        }
    }
    return 0;
}

/*
 * Declare system variables
 */
static realtype engine_time = 0;        /* Engine time */
static realtype engine_time_last = 0;   /* Previous engine time */
static realtype engine_pace = 0;
static realtype engine_realtime = 0;
static realtype engine_starttime = 0;
static realtype rootfinding_threshold = 0;
static long engine_evaluations = 0;
static long engine_steps = 0;

/*
 * Declare model and protocol variables
 */
Model model;                /* A model object */
ESys epacing;               /* Event-based pacing system */
FSys fpacing;               /* Fixed-form pacing system */

/*
 * Checks whether the model state matches that defined by the given arguments.
 *
 * Arguments
 *  t : Time
 *  p : Pacing variable
 *  y : The current state
 *  user_data : Model parameters (initial states are ignored)
 *
 * Returns 1 if the model is in the given state, 0 otherwise.
 */
static int
in_state(realtype t, realtype p, N_Vector y, void *user_data)
{
    /*TODO*/

    return 0;
}


/*
 * Right-hand-side function of the model ODE
 *
 *  realtype t      Current time
 *  N_Vector y      The current state values
 *  N_Vector ydot   Space to store the calculated derivatives in
 *  void* user_data Extra data (contains the sensitivity parameter values)
 *
 */
static int
rhs(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    FSys_Flag flag_fpacing;
    UserData fdata;
    int i;

    /* Recast user data */
    fdata = (UserData) user_data;

    /* Fixed-form pacing? Then look-up correct value of pacing variable! */
    if (fpacing != NULL) {
        engine_pace = FSys_GetLevel(fpacing, t, &flag_fpacing);
        if (flag_fpacing != FSys_OK) { /* This should never happen */
            FSys_SetPyErr(flag_fpacing);
            return -1;  /* Negative value signals irrecoverable error to CVODE */
        }
    }

    /*if (!in_state(t, engine_pace, y, user_data)) {*/
    if (1) {
        /* Update model state */

        /* Set time, pace, evaluations and realtime */
        engine_evaluations++;
        Model_UpdateBoundVariables(
            model, t, engine_pace, engine_realtime, engine_evaluations);

        /* Unpack parameters */
        if (model->n_parameters) {
            for (i=0; i<model->n_parameters; i++) {
                model->parameters[i] = fdata->p[i];
            }
            Model_UpdateParameterDerivedVariables(model);
        }

        /* Unpack states */
        for (i=0; i<model->n_states; i++) {
            model->states[i] = NV_Ith_S(y, i);
        }

        /* Calculate state derivatives */
        Model_UpdateDerivatives(model);
    }

    /* Fill ydot and return */
    for (i=0; i<model->n_states; i++) {
        NV_Ith_S(ydot, i) = model->derivatives[i];
    }

    return 0;
}

/*
 * Calculate sensitivities for output (starting from sensitivities of states)
 *
 * Assumes that all intermediary variables, derivatives, and sensitivities of
 * the derivatives have already been calculated.
 *
 *  realtype t      Current time
 *  N_Vector y      The current state values
 *  N_Vector ydot   The current state value derivatives
 *  N_Vector* yS    The current values of the sensitivity vectors (see below)
 *  N_Vector* ySdot Space to store any calculated sensitivities of derivatives
 *  void* user_data Extra data (contains the sensitivity parameter values)
 */
static int
calculate_sensitivity_outputs(realtype t, N_Vector y, N_Vector ydot,
                              N_Vector* yS, N_Vector* ySdot, void* user_data)
{
    int i, j, k;

    /* TODO: Memoisation */
    rhs(t, y, ydot, user_data);

    /* Unpack state sensitivities */
    for (i=0; i<model->n_sens; i++) {
        for (j=0; j<model->n_states; j++) {
            model->s_states[i * model->n_states + j] = NV_Ith_S(yS[i], j) * 0;
        }
    }

    /* Sensitivity w.r.t. ina.gNa */
    model->S0_V_INa = (((1.0 * pow(model->Y_m, 3.0) + model->P_gNa * (3.0 * pow(model->Y_m, 2.0) * model->S0_Y_m)) * model->Y_h + model->P_gNa * pow(model->Y_m, 3.0) * model->S0_Y_h) * model->Y_j + model->P_gNa * pow(model->Y_m, 3.0) * model->Y_h * model->S0_Y_j) * (model->Y_V - model->C_ENa) + model->P_gNa * pow(model->Y_m, 3.0) * model->Y_h * model->Y_j * model->S0_Y_V;

    /* Sensitivity w.r.t. ica.gCa */
    model->S1_V_INa = ((model->P_gNa * (3.0 * pow(model->Y_m, 2.0) * model->S1_Y_m) * model->Y_h + model->P_gNa * pow(model->Y_m, 3.0) * model->S1_Y_h) * model->Y_j + model->P_gNa * pow(model->Y_m, 3.0) * model->Y_h * model->S1_Y_j) * (model->Y_V - model->C_ENa) + model->P_gNa * pow(model->Y_m, 3.0) * model->Y_h * model->Y_j * model->S1_Y_V;

    /* Sensitivity w.r.t. init(ica.Ca_i) */
    model->S2_V_INa = ((model->P_gNa * (3.0 * pow(model->Y_m, 2.0) * model->S2_Y_m) * model->Y_h + model->P_gNa * pow(model->Y_m, 3.0) * model->S2_Y_h) * model->Y_j + model->P_gNa * pow(model->Y_m, 3.0) * model->Y_h * model->S2_Y_j) * (model->Y_V - model->C_ENa) + model->P_gNa * pow(model->Y_m, 3.0) * model->Y_h * model->Y_j * model->S2_Y_V;


    /* Fill ySdot and return */
    for (i=0; i<model->n_sens; i++) {
        for (j=0; j<model->n_states; j++) {
            NV_Ith_S(ySdot[i], j) = model->s_derivatives[i * model->n_states + j];
        }
    }

    return 0;
}

/*
 * Root finding function
 */
static int
root_finding(realtype t, N_Vector y, realtype *gout, void *user_data)
{
    gout[0] = NV_Ith_S(y, 0) - rootfinding_threshold;
    return 0;
}

/*
 * Settings
 */
static double abs_tol = 1e-6; /* The absolute tolerance */
static double rel_tol = 1e-4; /* The relative tolerance */
static double dt_max = 0;     /* The maximum step size (0.0 for none) */
static double dt_min = 0;     /* The minimum step size (0.0 for none) */

/*
 * Change the tolerance settings
 */
static PyObject*
sim_set_tolerance(PyObject *self, PyObject *args)
{
    /* Check input arguments */
    double tabs, trel;
    if (!PyArg_ParseTuple(args, "dd", &tabs, &trel)) {
        PyErr_SetString(PyExc_Exception, "Expected input arguments: abs_tol(float), rel_tol(float).");
        return 0;
    }
    abs_tol = tabs;
    rel_tol = trel;
    Py_RETURN_NONE;
}

/*
 * Change the maximum step size (0 for none)
 */
static PyObject*
sim_set_max_step_size(PyObject *self, PyObject *args)
{
    /* Check input arguments */
    double tmax;
    if (!PyArg_ParseTuple(args, "d", &tmax)) {
        PyErr_SetString(PyExc_Exception, "Expected input argument: tmax(float).");
        return 0;
    }
    dt_max = tmax;
    Py_RETURN_NONE;
}

/*
 * Change the minimum step size (0 for none)
 */
static PyObject*
sim_set_min_step_size(PyObject *self, PyObject *args)
{
    /* Check input arguments */
    double tmin;
    if (!PyArg_ParseTuple(args, "d", &tmin)) {
        PyErr_SetString(PyExc_Exception, "Expected input argument: tmin(float).");
        return 0;
    }
    dt_min = tmin;
    Py_RETURN_NONE;
}

/*
 * Add a variable to the logging lists. Returns 1 if successful
 */
static int
log_add(PyObject* log_dict, PyObject** logs, realtype** vars, int i, const char* name, const realtype* var)
{
    /* See first use of log_add for notes on unicode */
    int added = 0;
    PyObject* key = PyUnicode_FromString(name);
    if (PyDict_Contains(log_dict, key)) {
        logs[i] = PyDict_GetItem(log_dict, key);
        vars[i] = (realtype*)var;
        added = 1;
    }
    Py_DECREF(key);
    return added;
}

/*
 * Simulation variables
 */

int running = 0;        /* Running yes or no */

/* Input arguments */
double tmin;            /* The initial simulation time */
double tmax;            /* The final simulation time */
PyObject* state_in;     /* The initial state */
PyObject* state_out;    /* The final state */
PyObject* inputs;       /* A vector used to return the bound variables final values */
PyObject* literals;     /* A list of literal constant values */
PyObject* parameters;   /* A list of parameter values */
PyObject* eprotocol;    /* An event-based pacing protocol */
PyObject* fprotocol;    /* A fixed-form pacing protocol */
PyObject* log_dict;     /* The log dict */
double log_interval;    /* Periodic logging: The log interval (0 to disable) */
PyObject* log_times;    /* Point-list logging: List of points (None to disable) */
PyObject* sens_list;    /* List to store sensitivities in (or None if not enabled) */
PyObject* root_list;    /* List to store found roots in (or None if not enabled) */
double root_threshold;  /* Threshold to use for root finding */
PyObject* benchtime;    /* Callable time() function or None */

/* Next simulation halting point */
double tnext;

/* Checking for repeated zero size steps */
int zero_step_count;
int max_zero_step_count = 500;   /* Increased this from 50 */

/* CVode objects */
void *cvode_mem;     /* The memory used by the solver */
N_Vector y;          /* Stores the current position y */
N_Vector y_log;      /* Used to store y when logging */
N_Vector dy_log;     /* Used to store dy when logging */
N_Vector* sy;        /* Vector of state sensitivities, 1 per variable in `variables` */
N_Vector* sy_log;    /* Used to store y sensitivities when logging */
N_Vector* sdy_log;   /* Used to store dy sensitivities when logging */
N_Vector y_last;     /* Used to store previous value of y for error handling */
UserData udata;      /* UserData struct, used to pass in parameters */

#if MYOKIT_SUNDIALS_VERSION >= 30000
SUNMatrix sundense_matrix;          /* Dense matrix for linear solves */
SUNLinearSolver sundense_solver;    /* Linear solver object */
#endif

/* Parameter/initial-condition scales */
realtype* pbar;             /* One number per parameter/initial condition, giving something in the expected magnitude of the param/init */

/* Root finding */
int* rootsfound;     /* Used to store found roots */

/* Logging */
PyObject** logs;            /* An array of pointers to a PyObject */
realtype** vars;            /* An array of pointers to realtype */
int n_vars;                 /* Number of logging variables */
int log_bound;              /* True if logging bound variables */
int log_inter;              /* True if logging intermediary variables */
int log_deriv;              /* True if logging derivatives */
PyObject* list_update_str;  /* PyUnicode, used to call "append" method */
Py_ssize_t ilog;            /* Periodic/point-list logging: Index of next point */
double tlog;                /* Periodic/point-list logging: Next point */
int dynamic_logging;        /* True if logging every point. */

/*
 * Cleans up after a simulation
 */
static PyObject*
sim_clean()
{
    if (running != 0) {
        /* Logging */
        free(vars); vars = NULL;
        free(logs); logs = NULL;
        if (list_update_str != NULL) {
            Py_XDECREF(list_update_str); list_update_str = NULL;
        }

        /* Root finding */
        free(rootsfound); rootsfound = NULL;

        /* CVode arrays */
        if (y != NULL) { N_VDestroy_Serial(y); y = NULL; }
        if (dy_log != NULL) { N_VDestroy_Serial(dy_log); dy_log = NULL; }
        if (model != NULL && model->is_ode && !dynamic_logging) {
            if (y_log != NULL) { N_VDestroy_Serial(y_log); y_log = NULL; }
        }
        if (model != NULL && model->n_sens) {
            if (sy != NULL) { N_VDestroyVectorArray(sy, model->n_sens); sy = NULL; }
            if (sdy_log != NULL) { N_VDestroyVectorArray(sdy_log, model->n_sens); sdy_log = NULL; }
            if (model->is_ode && !dynamic_logging) {
                if (sy_log != NULL) { N_VDestroyVectorArray(sy_log, model->n_sens); sy_log = NULL; }
            }
        }

        /* CVode objects */
        CVodeFree(&cvode_mem); cvode_mem = NULL;
        #if MYOKIT_SUNDIALS_VERSION >= 30000
        SUNLinSolFree(sundense_solver); sundense_solver = NULL;
        SUNMatDestroy(sundense_matrix); sundense_matrix = NULL;
        #endif

        /* Free user data and parameter scale array*/
        free(pbar);
        if (udata != NULL) {
            free(udata->p);
            free(udata); udata = NULL;
        }

        /* Free pacing system space */
        ESys_Destroy(epacing); epacing = NULL;
        FSys_Destroy(fpacing); fpacing = NULL;

        /* Free model space */
        Model_Destroy(model); model = NULL;

        /* No longer running */
        running = 0;
    }

    /* Return 0, allowing the construct
        PyErr_SetString(PyExc_Exception, "Oh noes!");
        return sim_clean()
       to terminate a python function. */
    return 0;
}
static PyObject*
py_sim_clean(PyObject *self, PyObject *args)
{
    sim_clean();
    Py_RETURN_NONE;
}

/*
 * Initialise a run.
 * Called by the Python code's run(), followed by several calls to sim_step().
 */
static PyObject*
sim_init(PyObject *self, PyObject *args)
{
    int flag;
    int flag_cvode;
    Model_Flag flag_model;
    ESys_Flag flag_epacing;
    FSys_Flag flag_fpacing;

    int i, j;
    int log_first_point;

    Py_ssize_t pos;
    PyObject *flt;
    PyObject *key;
    PyObject* ret;
    PyObject *value;
    PyObject *l1, *l2;

    #ifndef SUNDIALS_DOUBLE_PRECISION
    PyErr_SetString(PyExc_Exception, "Sundials must be compiled with double precision.");
    /* No memory freeing is needed here, return directly */
    return 0;
    #endif

    /* Check if already running */
    if (running) {
        PyErr_SetString(PyExc_Exception, "Simulation already initialized.");
        return 0;
    }

    /* Set all pointers used in sim_clean to null */
    /* Model and pacing */
    model = NULL;
    epacing = NULL;
    fpacing = NULL;
    /* User data and parameter scaling */
    udata = NULL;
    pbar = NULL;
    /* CVode arrays */
    y = NULL;
    y_log = NULL;
    dy_log = NULL;
    sy = NULL;
    sy_log = NULL;
    sdy_log = NULL;
    /* Logging */
    logs = NULL;
    vars = NULL;
    list_update_str = NULL;
    log_times = NULL;
    /* Root finding */
    rootsfound = NULL;
    /* CVode objects */
    cvode_mem = NULL;
    #if MYOKIT_SUNDIALS_VERSION >= 30000
        sundense_matrix = NULL;
        sundense_solver = NULL;
    #endif

    /* Check input arguments     0123456789ABCDEF*/
    if (!PyArg_ParseTuple(args, "ddOOOOOOOOdOOOdO",
            &tmin,              /* 0. Float: initial time */
            &tmax,              /* 1. Float: final time */
            &state_in,          /* 2. List: initial state */
            &state_out,         /* 3. List: store final state here */
            &inputs,            /* 4. List: store final bound variables here */
            &literals,          /* 5. List: literal constant values */
            &parameters,        /* 6. List: parameter values */
            &eprotocol,         /* 7. Event-based protocol */
            &fprotocol,         /* 8. Fixed-form protocol (tuple) */
            &log_dict,          /* 9. DataLog */
            &log_interval,      /* A. Float: log interval, or 0 */
            &log_times,         /* B. List of logging times, or None */
            &sens_list,         /* C. List to store sensitivities in */
            &root_list,         /* D. List to store roots in or None */
            &root_threshold,    /* E. Float: root-finding threshold */
            &benchtime          /* F. Callable to obtain system time */
            )) {
        PyErr_SetString(PyExc_Exception, "Incorrect input arguments.");
        /* Nothing allocated yet, no pyobjects _created_, return directly */
        return 0;
    }

    /* Now officialy running :) */
    running = 1;

    /*************************************************************************
    From this point on, no more direct returning! Use sim_clean()

    To check if this list is still up to date manually search for cvode
    and python stuff. To find what to free() search for "alloc("
    Initialize all to NULL so that free() will work without errors.

    Notes:
    1. Functions like PyList_New and PyDict_New create a new object with a
       refcount of 1. They pass on the ownership of this reference to the
       caller, IE they return the reference and it becomes the caller's
       responsibility to call PyDECREF
    2. Functions like PyList_Append and PyDict_SetItem create a new reference
       to the items you pass them, IE they increase the ref count and will
       decrease it when they're done with it. This means that you retain
       ownership of your own reference to this items and will also need to
       call decref when you're done with them.
    3. PyList_SetItem and PyTuple_SetItem are exceptions to the rule: they
       "steal" a reference to the item you pass into them. This means they do
       not increase the refcount of the item, but _do_ decrease it when they
       themselves are destructed.
       This _only_ holds for the SetItem functions, and _only_ for list and
       tuple.
       The reasonining behind this is that it's a very common scenario for
       populating lists and tuples.
    4. PyList_GetItem and PyTuple_GetItem are exceptions to the rule: they
       return a "borrowed" reference to an item. This means you should never
       decref them!
       This _only_ holds for list and tuple.
    5. When you return a newly created reference from a function, you pass on
       the ownership of that reference to the calling function. This means you
       don't have to call DECREF on the return value of a function.
    6. References passed _into_ your function as arguments are _borrowed_:
       Their refcount doesn't change and you don't have to increase or decrease
       it. The object they point to is guaranteed to exist for as long as your
       function runs.

    Result:
    A. The log and protocol objects passed to this function are borrowed
       references: no need to change the reference count.
    B. The PyFloat objects that are created have refcount 1. They're added to
       the lists using append, which increases their refcount. So they should
       be decref'd after appending.
    C. The time float that is created has refcount 1. It's ownership is passed
       on to the calling function. No need to decref.
    D. The PyFloat objects in this list are added using PyList_SetItem which
       steals ownership: No need to decref.
    */

    /* Create model */
    model = Model_Create(&flag_model);
    if (flag_model != Model_OK) { Model_SetPyErr(flag_model); return sim_clean(); }

    /* Create state vector */
    y = N_VNew_Serial(model->n_states);
    if (check_cvode_flag((void*)y, "N_VNew_Serial", 0)) {
        PyErr_SetString(PyExc_Exception, "Failed to create state vector.");
        return sim_clean();
    }

    /* Create state vector copy for error handling */
    y_last = N_VNew_Serial(model->n_states);
    if (check_cvode_flag((void*)y_last, "N_VNew_Serial", 0)) {
        PyErr_SetString(PyExc_Exception, "Failed to create last-state vector.");
        return sim_clean();
    }

    /* Create sensitivity vector array */
    if (model->n_sens) {
        sy = N_VCloneVectorArray(model->n_sens, y);
        if (check_cvode_flag((void*)sy, "N_VCloneVectorArray", 0)) {
            PyErr_SetString(PyExc_Exception, "Failed to allocate space to store sensitivities.");
            return sim_clean();
        }
    }

   /* Determine if dynamic logging is being used (or if it's periodic/point-list logging) */
    dynamic_logging = (log_interval <= 0 && log_times == Py_None);

    /* Create state vector for logging.
       If we log every visited point (dynamic logging), we can simply log
       values from y as used for the integration. But if we need to interpolate
       times in the past, we need a vector y (and sy if using sensitivities) to
       store the state at time of logging in. */
    if (dynamic_logging || !model->is_ode) {
        /* Dynamic logging or cvode-free mode: don't interpolate,
           so let y_log point to y */
        y_log = y;
        sy_log = sy;
    } else {
        /* Logging at fixed points:
           Keep y_log as a separate N_Vector for cvode interpolation */
        y_log = N_VNew_Serial(model->n_states);
        if (check_cvode_flag((void*)y_log, "N_VNew_Serial", 0)) {
            PyErr_SetString(PyExc_Exception, "Failed to create state vector for logging.");
            return sim_clean();
        }
        if (model->n_sens) {
            sy_log = N_VCloneVectorArray(model->n_sens, y);
            if (check_cvode_flag((void*)sy_log, "N_VCloneVectorArray", 0)) {
                PyErr_SetString(PyExc_Exception, "Failed to create state sensitivity vector array for logging.");
                return sim_clean();
            }
        }
    }

    /* Create derivative vector for logging.
       In both logging modes, the need arises to call rhs() manually, and so we
       need to create a dy_log vector to pass in (and an sdy_log vector array).
    */
    dy_log = N_VNew_Serial(model->n_states);
    if (check_cvode_flag((void*)dy_log, "N_VNew_Serial", 0)) {
        PyErr_SetString(PyExc_Exception, "Failed to create derivatives vector for logging.");
        return sim_clean();
    }
    if (model->n_sens) {
        sdy_log = N_VCloneVectorArray(model->n_sens, y);
        if (check_cvode_flag((void*)sdy_log, "N_VCloneVectorArray", 0)) {
            PyErr_SetString(PyExc_Exception, "Failed to create derivatives sensitivity vector array for logging.");
            return sim_clean();
        }
    }

    /* Set values of literals */
    if (!PyList_Check(literals)) {
        PyErr_SetString(PyExc_Exception, "'literals' must be a list.");
        return sim_clean();
    }
    for (i=0; i<model->n_literals; i++) {
        flt = PyList_GetItem(literals, i);    /* Don't decref */
        if (!PyFloat_Check(flt)) {
            char errstr[200];
            sprintf(errstr, "Item %d in literal vector is not a float.", i);
            PyErr_SetString(PyExc_Exception, errstr);
            return sim_clean();
        }
        model->literals[i] = PyFloat_AsDouble(flt);
    }

    /* Set calculated constants */
    Model_UpdateLiteralDerivedVariables(model);

    /* The model's parameter values will be set when rhs is called. */

    /* Set initial state values */
    if (!PyList_Check(state_in)) {
        PyErr_SetString(PyExc_Exception, "'state_in' must be a list.");
        return sim_clean();
    }
    for (i=0; i<model->n_states; i++) {
        flt = PyList_GetItem(state_in, i);    /* Don't decref! */
        if (!PyFloat_Check(flt)) {
            char errstr[200];
            sprintf(errstr, "Item %d in state vector is not a float.", i);
            PyErr_SetString(PyExc_Exception, errstr);
            return sim_clean();
        }
        NV_Ith_S(y, i) = PyFloat_AsDouble(flt);
        NV_Ith_S(y_last, i) = NV_Ith_S(y, i);
        /*if (!dynamic_logging) {
            NV_Ith_S(y_log, i) = NV_Ith_S(y, i);
        }*/
    }

    /* Set initial sensitivities to zero, or to 1 for initial conditions */
    /*TODO: GET THESE FROM THE PYTHON CODE */
    for (i=0; i<model->n_sens; i++) {
        N_VConst(RCONST(0.0), sy[i]);
    }
    for (i=0; i<model->n_initials; i++) {
        NV_Ith_S(sy[model->n_parameters + i], model->i_initials[i]) = 1.0;
    }

    /* Create UserData and vector to hold parameter and initial condition values for cvode */
    udata = (UserData)malloc(sizeof *udata);
    if (udata == 0) {
        PyErr_SetString(PyExc_Exception, "Unable to create user data object to store parameter values.");
        return sim_clean();
    }
    udata->p = (realtype*)malloc(sizeof(realtype) * model->n_sens);
    if (udata->p == 0) {
        PyErr_SetString(PyExc_Exception, "Unable to allocate space to store parameter values.");
        return sim_clean();
    }

    /* Add parameter values to UserData */
    if (!PyList_Check(parameters)) {
        PyErr_SetString(PyExc_Exception, "'parameters' must be a list.");
        return sim_clean();
    }
    for (i=0; i<model->n_parameters; i++) {
        flt = PyList_GetItem(parameters, i);    /* Don't decref */
        if (!PyFloat_Check(flt)) {
            char errstr[200];
            sprintf(errstr, "Item %d in parameter vector is not a float.", i);
            PyErr_SetString(PyExc_Exception, errstr);
            return sim_clean();
        }
        udata->p[i] = PyFloat_AsDouble(flt);
    }

    /* Add initial values to UserData */
    for (i=0; i<model->n_initials; i++) {
        udata->p[model->n_parameters + i] = NV_Ith_S(y, model->i_initials[i]);
    }

    /* Create parameter scaling vector, for error control */
    pbar = (realtype*)malloc(sizeof(realtype) * model->n_sens);
    if (pbar == 0) {
        PyErr_SetString(PyExc_Exception, "Unable to allocate space to store parameter scales.");
        return sim_clean();
    }
    for (i=0; i<model->n_sens; i++) {
        pbar[i] = (udata->p[i] == 0.0 ? 1.0 : fabs(udata->p[i]));
    }

    /* Root finding list of integers (only contains 1 int...) */
    rootsfound = (int*)malloc(sizeof(int)*1);

    /* Reset evaluation count */
    engine_evaluations = 0;

    /* Reset step count */
    engine_steps = 0;

    /* Zero step tracking */
    zero_step_count = 0;

    /* Check output list */
    if (!PyList_Check(state_out)) {
        PyErr_SetString(PyExc_Exception, "'state_out' must be a list.");
        return sim_clean();
    }

    /* Check for loss-of-precision issue in periodic logging */
    if (log_interval > 0) {
        if (tmax + log_interval == tmax) {
            PyErr_SetString(PyExc_Exception, "Log interval is too small compared to tmax; issue with numerical precision: float(tmax + log_interval) = float(tmax).");
            return sim_clean();
        }
    }

    /* Set up logging */
    log_inter = 0;
    log_bound = 0;
    n_vars = PyDict_Size(log_dict);
    logs = (PyObject**)malloc(sizeof(PyObject*)*n_vars);
    vars = (realtype**)malloc(sizeof(realtype*)*n_vars);
    i = 0;

    /* Note: The variable names are all ascii compatible
       In Python2, they are stored in logs as either unicode or bytes
       In Python3, they are stored exclusively as unicode
       However, in Python2 b'name' matches u'name' so this is ok (in Python it
        does not, but we always use unicode so it's ok).
       The strategy here will be to convert these C-strings to unicode
        inside log_add, before the comparison. */

    /* Check states */
    i += log_add(log_dict, logs, vars, i, "membrane.V", &NV_Ith_S(y_log, 0));
    i += log_add(log_dict, logs, vars, i, "ina.m", &NV_Ith_S(y_log, 1));
    i += log_add(log_dict, logs, vars, i, "ina.h", &NV_Ith_S(y_log, 2));
    i += log_add(log_dict, logs, vars, i, "ina.j", &NV_Ith_S(y_log, 3));
    i += log_add(log_dict, logs, vars, i, "ica.d", &NV_Ith_S(y_log, 4));
    i += log_add(log_dict, logs, vars, i, "ica.f", &NV_Ith_S(y_log, 5));
    i += log_add(log_dict, logs, vars, i, "ik.x", &NV_Ith_S(y_log, 6));
    i += log_add(log_dict, logs, vars, i, "ica.Ca_i", &NV_Ith_S(y_log, 7));


    /* Check derivatives */
    j = i;
    i += log_add(log_dict, logs, vars, i, "dot(membrane.V)", &NV_Ith_S(dy_log, 0));
    i += log_add(log_dict, logs, vars, i, "dot(ina.m)", &NV_Ith_S(dy_log, 1));
    i += log_add(log_dict, logs, vars, i, "dot(ina.h)", &NV_Ith_S(dy_log, 2));
    i += log_add(log_dict, logs, vars, i, "dot(ina.j)", &NV_Ith_S(dy_log, 3));
    i += log_add(log_dict, logs, vars, i, "dot(ica.d)", &NV_Ith_S(dy_log, 4));
    i += log_add(log_dict, logs, vars, i, "dot(ica.f)", &NV_Ith_S(dy_log, 5));
    i += log_add(log_dict, logs, vars, i, "dot(ik.x)", &NV_Ith_S(dy_log, 6));
    i += log_add(log_dict, logs, vars, i, "dot(ica.Ca_i)", &NV_Ith_S(dy_log, 7));

    log_deriv = (i > j);

    /* Check bound variables */
    j = i;
    i += log_add(log_dict, logs, vars, i, "engine.time", &model->V_time);
    i += log_add(log_dict, logs, vars, i, "engine.pace", &model->V_pace);

    log_bound = (i > j);

    /* Remaining variables will require an extra rhs() call to evaluate their
       values at every log point */
    j = i;
    i += log_add(log_dict, logs, vars, i, "membrane.i_ion", &model->V_i_ion);
    i += log_add(log_dict, logs, vars, i, "membrane.i_stim", &model->V_i_stim);
    i += log_add(log_dict, logs, vars, i, "ik.x.alpha", &model->V_ik_x_alpha);
    i += log_add(log_dict, logs, vars, i, "ik.x.beta", &model->V_ik_x_beta);
    i += log_add(log_dict, logs, vars, i, "ik.xi", &model->V_xi);
    i += log_add(log_dict, logs, vars, i, "ik.IK", &model->V_IK);
    i += log_add(log_dict, logs, vars, i, "ina.a", &model->V_a);
    i += log_add(log_dict, logs, vars, i, "ina.m.alpha", &model->V_ina_m_alpha);
    i += log_add(log_dict, logs, vars, i, "ina.m.beta", &model->V_ina_m_beta);
    i += log_add(log_dict, logs, vars, i, "ina.h.alpha", &model->V_ina_h_alpha);
    i += log_add(log_dict, logs, vars, i, "ina.h.beta", &model->V_ina_h_beta);
    i += log_add(log_dict, logs, vars, i, "ina.j.alpha", &model->V_ina_j_alpha);
    i += log_add(log_dict, logs, vars, i, "ina.j.beta", &model->V_ina_j_beta);
    i += log_add(log_dict, logs, vars, i, "ina.INa", &model->V_INa);
    i += log_add(log_dict, logs, vars, i, "ikp.Kp", &model->V_Kp);
    i += log_add(log_dict, logs, vars, i, "ikp.IKp", &model->V_IKp);
    i += log_add(log_dict, logs, vars, i, "ica.E", &model->V_ica_E);
    i += log_add(log_dict, logs, vars, i, "ica.d.alpha", &model->V_ica_d_alpha);
    i += log_add(log_dict, logs, vars, i, "ica.d.beta", &model->V_ica_d_beta);
    i += log_add(log_dict, logs, vars, i, "ica.f.alpha", &model->V_ica_f_alpha);
    i += log_add(log_dict, logs, vars, i, "ica.f.beta", &model->V_ica_f_beta);
    i += log_add(log_dict, logs, vars, i, "ica.ICa", &model->V_ICa);
    i += log_add(log_dict, logs, vars, i, "ik1.g", &model->V_g);
    i += log_add(log_dict, logs, vars, i, "ik1.g.alpha", &model->V_ik1_g_alpha);
    i += log_add(log_dict, logs, vars, i, "ik1.g.beta", &model->V_ik1_g_beta);
    i += log_add(log_dict, logs, vars, i, "ik1.IK1", &model->V_IK1);
    i += log_add(log_dict, logs, vars, i, "ib.Ib", &model->V_Ib);

    log_inter = (i > j);

    /* Check if log contained extra variables */
    if (i != n_vars) {
        PyErr_SetString(PyExc_Exception, "Unknown variables found in logging dictionary.");
        return sim_clean();
    }

    /* Check logging list for sensitivities */
    if (model->n_sens) {
        if (!PyList_Check(sens_list)) {
            PyErr_SetString(PyExc_Exception, "'sens_list' must be a list.");
            return sim_clean();
        }
    }

    /* Set up event-based pacing */
    if (eprotocol != Py_None) {
        epacing = ESys_Create(&flag_epacing);
        if (flag_epacing != ESys_OK) { ESys_SetPyErr(flag_epacing); return sim_clean(); }
        flag_epacing = ESys_Populate(epacing, eprotocol);
        if (flag_epacing != ESys_OK) { ESys_SetPyErr(flag_epacing); return sim_clean(); }
        flag_epacing = ESys_AdvanceTime(epacing, tmin);
        if (flag_epacing != ESys_OK) { ESys_SetPyErr(flag_epacing); return sim_clean(); }
        tnext = ESys_GetNextTime(epacing, &flag_epacing);
        engine_pace = ESys_GetLevel(epacing, &flag_epacing);
        tnext = (tnext < tmax) ? tnext : tmax;
    } else {
        tnext = tmax;
    }

    /* Set up fixed-form pacing */
    if (eprotocol == Py_None && fprotocol != Py_None) {
        /* Check 'protocol' is tuple (times, values) */
        if (!PyTuple_Check(fprotocol)) {
            PyErr_SetString(PyExc_Exception, "Fixed-form pacing protocol should be tuple or None.");
            return sim_clean();
        }
        if (PyTuple_Size(fprotocol) != 2) {
            PyErr_SetString(PyExc_Exception, "Fixed-form pacing protocol tuple should have size 2.");
            return sim_clean();
        }
        /* Create fixed-form pacing object and populate */
        fpacing = FSys_Create(&flag_fpacing);
        if (flag_fpacing != FSys_OK) { FSys_SetPyErr(flag_fpacing); return sim_clean(); }
        flag_fpacing = FSys_Populate(fpacing,
            PyTuple_GetItem(fprotocol, 0),  /* Borrowed, no decref */
            PyTuple_GetItem(fprotocol, 1));
        if (flag_fpacing != FSys_OK) { FSys_SetPyErr(flag_fpacing); return sim_clean(); }
    }

    /* Set simulation starting time */
    engine_time = tmin;

    /* Check dt_max and dt_min */
    if (dt_max < 0) dt_max = 0.0;
    if (dt_min < 0) dt_min = 0.0;

    /* Create solver
     * Using Backward differentiation and Newton iteration */
    if (model->is_ode) {
        #if MYOKIT_SUNDIALS_VERSION >= 40000
            cvode_mem = CVodeCreate(CV_BDF);
        #else
            cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
        #endif
        if (check_cvode_flag((void*)cvode_mem, "CVodeCreate", 0)) return sim_clean();

        /* Initialise solver memory, specify the rhs */
        flag_cvode = CVodeInit(cvode_mem, rhs, engine_time, y);
        if (check_cvode_flag(&flag_cvode, "CVodeInit", 1)) return sim_clean();

        /* Set absolute and relative tolerances */
        flag_cvode = CVodeSStolerances(cvode_mem, RCONST(rel_tol), RCONST(abs_tol));
        if (check_cvode_flag(&flag_cvode, "CVodeSStolerances", 1)) return sim_clean();

        /* Set a maximum step size (or 0.0 for none) */

        flag_cvode = CVodeSetMaxStep(cvode_mem, dt_max);
        if (check_cvode_flag(&flag_cvode, "CVodeSetmaxStep", 1)) return sim_clean();

        /* Set a minimum step size (or 0.0 for none) */
        flag_cvode = CVodeSetMinStep(cvode_mem, dt_min);
        if (check_cvode_flag(&flag_cvode, "CVodeSetminStep", 1)) return sim_clean();

        #if MYOKIT_SUNDIALS_VERSION >= 30000
            /* Create dense matrix for use in linear solves */
            sundense_matrix = SUNDenseMatrix(model->n_states, model->n_states);
            if (check_cvode_flag((void *)sundense_matrix, "SUNDenseMatrix", 0)) return sim_clean();

            /* Create dense linear solver object with matrix */
            sundense_solver = SUNDenseLinearSolver(y, sundense_matrix);
            if (check_cvode_flag((void *)sundense_solver, "SUNDenseLinearSolver", 0)) return sim_clean();

            /* Attach the matrix and solver to cvode */
            flag_cvode = CVDlsSetLinearSolver(cvode_mem, sundense_solver, sundense_matrix);
            if (check_cvode_flag(&flag_cvode, "CVDlsSetLinearSolver", 1)) return sim_clean();
        #else
            /* Create dense matrix for use in linear solves */
            flag_cvode = CVDense(cvode_mem, model->n_states);
            if (check_cvode_flag(&flag_cvode, "CVDense", 1)) return sim_clean();
        #endif

        /* Activate forward sensitivity computations */
        if (model->n_sens) {
            /* TODO: NULL here is the place to insert a user function to calculate the
               RHS of the sensitivity ODE */
            /*flag_cvode = CVodeSensInit(cvode_mem, model->n_sens, CV_SIMULTANEOUS, rhs1, sy);*/
            flag_cvode = CVodeSensInit(cvode_mem, model->n_sens, CV_SIMULTANEOUS, NULL, sy);
            if (check_cvode_flag(&flag_cvode, "CVodeSensInit", 1)) return sim_clean();

            /* Attach user data */
            flag_cvode = CVodeSetUserData(cvode_mem, udata);
            if (check_cvode_flag(&flag_cvode, "CVodeSetUserData", 1)) return sim_clean();

            /* Set parameter scales used in tolerances */
            flag_cvode = CVodeSetSensParams(cvode_mem, udata->p, pbar, NULL);
            if (check_cvode_flag(&flag_cvode, "CVodeSetSensParams", 1)) return sim_clean();

            /* Set sensitivity tolerances calculating method (using pbar) */
            flag_cvode = CVodeSensEEtolerances(cvode_mem);
            if (check_cvode_flag(&flag_cvode, "CVodeSensEEtolerances", 1)) return sim_clean();
        }

    } /* if model.CVODE */

    /* Set string for updating lists/arrays using Python interface. */
    list_update_str = PyUnicode_FromString("append");

    /* Benchmarking? Then set engine_realtime to 0.0 */
    if (benchtime != Py_None) {
        /* Store initial time as 0 */
        engine_realtime = 0.0;
        /* Tell sim_step to set engine_starttime */
        engine_starttime = -1;
    }

    /* Set logging points */
    if (log_interval > 0) {

        /* Periodic logging */
        ilog = 0;
        tlog = tmin;

    } else if (log_times != Py_None) {

        /* Point-list logging */

        /* Check the log_times list */
        if (!PyList_Check(log_times)) {
            PyErr_SetString(PyExc_Exception, "'log_times' must be a list.");
            return sim_clean();
        }

        /* Read next log point off the list */
        ilog = 0;
        tlog = engine_time - 1;
        while(ilog < PyList_Size(log_times) && tlog < engine_time) {
            flt = PyList_GetItem(log_times, ilog); /* Borrowed */
            if (!PyFloat_Check(flt)) {
                PyErr_SetString(PyExc_Exception, "Entries in 'log_times' must be floats.");
                return sim_clean();
            }
            tlog = PyFloat_AsDouble(flt);
            ilog++;
            flt = NULL;
        }

        /* No points beyond engine_time? Then don't log any future points. */
        if (tlog < engine_time) {
            tlog = tmax + 1;
        }

    } else {

        /*
         * Dynamic logging
         *
         * Log the first entry, but only if not appending to an existing log.
         * This prevents points from appearing twice when a simulation with
         * dynamic logging is stopped and started.
         */

        /*  Check if the log is empty */
        log_first_point = 1;
        pos = 0;
        if (PyDict_Next(log_dict, &pos, &key, &value)) {
            /* Items found in dict, randomly selected list now in "value" */
            /* Both key and value are borrowed references, no need to decref */
            log_first_point = (PyObject_Size(value) <= 0);
        }

        /* If so, log the first point! */
        if (log_first_point) {
            rhs(engine_time, y, dy_log, udata);
            /* At this point, we have y(t), inter(t) and dy(t) */
            /* We've also loaded time(t) and pace(t) */
            for (i=0; i<n_vars; i++) {
                flt = PyFloat_FromDouble(*vars[i]);
                ret = PyObject_CallMethodObjArgs(logs[i], list_update_str, flt, NULL);
                Py_DECREF(flt);
                Py_XDECREF(ret);
                if (ret == NULL) {
                    flt = NULL;
                    PyErr_SetString(PyExc_Exception, "Call to append() failed on logging list.");
                    return sim_clean();
                }
            }
            flt = NULL;
            ret = NULL;

            if (model->n_sens) {
                /* Calculate sensitivities to output */

                /* TODO: Call the rhs function that sets sdy_log */
                calculate_sensitivity_outputs(engine_time, y, dy_log, sy, sdy_log, udata);

                /* Write sensitivity matrix to log */
                l1 = PyTuple_New(2);
                if (l1 == NULL) return sim_clean();

                l2 = PyTuple_New(model->n_sens);
                if (l2 == NULL) return sim_clean();

                flt = PyFloat_FromDouble(model->S0_Y_V);
                if (flt == NULL) return sim_clean();
                flag = PyTuple_SetItem(l2, 0, flt); /* Steals reference to flt */
                if (flag < 0) return sim_clean();
                flt = PyFloat_FromDouble(model->S1_Y_V);
                if (flt == NULL) return sim_clean();
                flag = PyTuple_SetItem(l2, 1, flt); /* Steals reference to flt */
                if (flag < 0) return sim_clean();
                flt = PyFloat_FromDouble(model->S2_Y_V);
                if (flt == NULL) return sim_clean();
                flag = PyTuple_SetItem(l2, 2, flt); /* Steals reference to flt */
                if (flag < 0) return sim_clean();

                flag = PyTuple_SetItem(l1, 0, l2); /* Steals reference to l2 */
                l2 = NULL; flt = NULL;

                l2 = PyTuple_New(model->n_sens);
                if (l2 == NULL) return sim_clean();

                flt = PyFloat_FromDouble(model->S0_V_INa);
                if (flt == NULL) return sim_clean();
                flag = PyTuple_SetItem(l2, 0, flt); /* Steals reference to flt */
                if (flag < 0) return sim_clean();
                flt = PyFloat_FromDouble(model->S1_V_INa);
                if (flt == NULL) return sim_clean();
                flag = PyTuple_SetItem(l2, 1, flt); /* Steals reference to flt */
                if (flag < 0) return sim_clean();
                flt = PyFloat_FromDouble(model->S2_V_INa);
                if (flt == NULL) return sim_clean();
                flag = PyTuple_SetItem(l2, 2, flt); /* Steals reference to flt */
                if (flag < 0) return sim_clean();

                flag = PyTuple_SetItem(l1, 1, l2); /* Steals reference to l2 */
                l2 = NULL; flt = NULL;

                flag = PyList_Append(sens_list, l1);
                Py_XDECREF(l1); l1 = NULL;
                if (flag < 0) return sim_clean();

            }
        }
    }

    /* Root finding enabled? (cvode-mode only) */
    if (model->is_ode && PySequence_Check(root_list)) {
        /* Set threshold */
        rootfinding_threshold = root_threshold;
        /* Initialize root function with 1 component */
        flag_cvode = CVodeRootInit(cvode_mem, 1, root_finding);
        if (check_cvode_flag(&flag_cvode, "CVodeRootInit", 1)) return sim_clean();
    }

    /* Done! */
    Py_RETURN_NONE;
}

/*
 * Takes the next steps in a simulation run
 */
static PyObject*
sim_step(PyObject *self, PyObject *args)
{
    ESys_Flag flag_epacing;
    int i;
    int steps_taken = 0;    /* Number of integration steps taken in this call */
    int flag_cvode;         /* CVode flag */
    int flag_root;          /* Root finding flag */
    int flag_reinit = 0;    /* Set if CVODE needs to be reset during a simulation step */
    int flag;               /* General flag for python operations */
    PyObject *flt, *ret;
    PyObject *l1, *l2;

    /*
     * Benchmarking? Then make sure start time is set.
     * This is handled here instead of in sim_init so it only includes time
     * taken performing steps, not time initialising memory etc.
     */
    if (benchtime != Py_None && engine_starttime < 0) {
        flt = PyObject_CallFunction(benchtime, "");
        if (!PyFloat_Check(flt)) {
            Py_XDECREF(flt); flt = NULL;
            PyErr_SetString(PyExc_Exception, "Call to benchmark time function didn't return float.");
            return sim_clean();
        }
        engine_starttime = PyFloat_AsDouble(flt);
        Py_DECREF(flt); flt = NULL;
    }

    /* Go! */
    while(1) {

        /* Back-up current y (no allocation, this is fast) */
        for (i=0; i<model->n_states; i++) {
            NV_Ith_S(y_last, i) = NV_Ith_S(y, i);
        }

        /* Store engine time before step */
        engine_time_last = engine_time;

        if (model->is_ode) {

            /* Take a single ODE step */
            flag_cvode = CVode(cvode_mem, tnext, y, &engine_time, CV_ONE_STEP);

            /* Check for errors */
            if (check_cvode_flag(&flag_cvode, "CVode", 1)) {
                /* Something went wrong... Set outputs and return */
                for (i=0; i<model->n_states; i++) {
                    PyList_SetItem(state_out, i, PyFloat_FromDouble(NV_Ith_S(y_last, i)));
                    /* PyList_SetItem steals a reference: no need to decref the double! */
                }
                PyList_SetItem(inputs, 0, PyFloat_FromDouble(engine_time));
                PyList_SetItem(inputs, 1, PyFloat_FromDouble(engine_pace));
                PyList_SetItem(inputs, 2, PyFloat_FromDouble(engine_realtime));
                PyList_SetItem(inputs, 3, PyFloat_FromDouble(engine_evaluations));
                return sim_clean();
            }

        } else {

            /* Just jump to next event */
            /* Note 1: To stay compatible with cvode-mode, don't jump to the
               next log time (if tlog < tnext) */
            /* Note 2: tnext can be infinity, so don't always jump there. */
            engine_time = (tmax > tnext) ? tnext : tmax;
            flag_cvode = CV_SUCCESS;

        }

        /* Check if progress is being made */
        if (engine_time == engine_time_last) {
            if (++zero_step_count >= max_zero_step_count) {
                char errstr[200];
                sprintf(errstr, "ZERO_STEP %f", engine_time);
                PyErr_SetString(PyExc_Exception, errstr);
                return sim_clean();
            }
        } else {
            /* Only count consecutive zero steps! */
            zero_step_count = 0;
        }

        /* Update step count */
        engine_steps++;

        /* If we got to this point without errors... */
        if ((flag_cvode == CV_SUCCESS) || (flag_cvode == CV_ROOT_RETURN)) {

            /* Interpolation and root finding */
            if (model->is_ode) {

                /* Next event time exceeded? */
                if (engine_time > tnext) {

                    /* Go back to engine_time=tnext */
                    flag_cvode = CVodeGetDky(cvode_mem, tnext, 0, y);
                    if (check_cvode_flag(&flag_cvode, "CVodeGetDky", 1)) return sim_clean();
                    if (model->n_sens) {
                        flag_cvode = CVodeGetSensDky(cvode_mem, tnext, 0, sy);
                        if (check_cvode_flag(&flag_cvode, "CVodeGetSensDky", 1)) return sim_clean();
                    }
                    engine_time = tnext;
                    /* Require reinit (after logging) */
                    flag_reinit = 1;

                } else {
                    /* Get current sensitivity vector */
                    if (model->n_sens) {
                        flag_cvode = CVodeGetSens(cvode_mem, &engine_time, sy);
                        if (check_cvode_flag(&flag_cvode, "CVodeGetSens", 1)) return sim_clean();
                    }

                    /* Root found */
                    if (flag_cvode == CV_ROOT_RETURN) {

                        /* Store found roots */
                        flag_root = CVodeGetRootInfo(cvode_mem, rootsfound);
                        if (check_cvode_flag(&flag_root, "CVodeGetRootInfo", 1)) return sim_clean();
                        flt = PyTuple_New(2);
                        PyTuple_SetItem(flt, 0, PyFloat_FromDouble(engine_time)); /* Steals reference, so this is ok */
                        PyTuple_SetItem(flt, 1, PyLong_FromLong(rootsfound[0]));
                        ret = PyObject_CallMethodObjArgs(root_list, list_update_str, flt, NULL);
                        Py_DECREF(flt); flt = NULL;
                        Py_XDECREF(ret);
                        if (ret == NULL) {
                            PyErr_SetString(PyExc_Exception, "Call to append() failed on root finding list.");
                            return sim_clean();
                        }
                        ret = NULL;
                    }
                }
            }

            /* Periodic logging or point-list logging */
            if (!dynamic_logging && engine_time > tlog) {
                /* Note: For periodic logging, the condition should be
                   `time > tlog` so that we log half-open intervals (i.e. the
                   final point should never be included). */

                /* Benchmarking? Then set engine_realtime */
                if (benchtime != Py_None) {
                    flt = PyObject_CallFunction(benchtime, "");
                    if (!PyFloat_Check(flt)) {
                        Py_XDECREF(flt); flt = NULL;
                        PyErr_SetString(PyExc_Exception, "Call to benchmark time function didn't return float.");
                        return sim_clean();
                    }
                    engine_realtime = PyFloat_AsDouble(flt) - engine_starttime;
                    Py_DECREF(flt); flt = NULL;
                }

                /* Log points */
                while (engine_time > tlog) {

                    /* Get interpolated y(tlog) */
                    if (model->is_ode) {
                        flag_cvode = CVodeGetDky(cvode_mem, tlog, 0, y_log);
                        if (check_cvode_flag(&flag_cvode, "CVodeGetDky", 1)) return sim_clean();
                        if (model->n_sens) {
                            flag_cvode = CVodeGetSensDky(cvode_mem, tlog, 0, sy_log);
                            if (check_cvode_flag(&flag_cvode, "CVodeGetSensDky", 1)) return sim_clean();
                        }
                    }
                    /* If cvode-free mode, the state can't change so we don't
                       need to do anything here */

                    /* Calculate intermediate variables & derivatives */
                    rhs(tlog, y_log, dy_log, udata);

                    /* Write to log */
                    for (i=0; i<n_vars; i++) {
                        flt = PyFloat_FromDouble(*vars[i]);
                        ret = PyObject_CallMethodObjArgs(logs[i], list_update_str, flt, NULL);
                        Py_DECREF(flt);
                        Py_XDECREF(ret);
                        if (ret == NULL) {
                            flt = NULL;
                            PyErr_SetString(PyExc_Exception, "Call to append() failed on logging list.");
                            return sim_clean();
                        }
                    }
                    ret = flt = NULL;

                    if (model->n_sens) {
                        /* Calculate sensitivities to output */
                        calculate_sensitivity_outputs(engine_time, y_log, dy_log, sy_log, sdy_log, udata);

                        /* Write sensitivity matrix to log */
                        l1 = PyTuple_New(2);
                        if (l1 == NULL) return sim_clean();

                        l2 = PyTuple_New(model->n_sens);
                        if (l2 == NULL) return sim_clean();

                        flt = PyFloat_FromDouble(model->S0_Y_V);
                        if (flt == NULL) return sim_clean();
                        flag = PyTuple_SetItem(l2, 0, flt); /* Steals reference to flt */
                        if (flag < 0) return sim_clean();
                        flt = PyFloat_FromDouble(model->S1_Y_V);
                        if (flt == NULL) return sim_clean();
                        flag = PyTuple_SetItem(l2, 1, flt); /* Steals reference to flt */
                        if (flag < 0) return sim_clean();
                        flt = PyFloat_FromDouble(model->S2_Y_V);
                        if (flt == NULL) return sim_clean();
                        flag = PyTuple_SetItem(l2, 2, flt); /* Steals reference to flt */
                        if (flag < 0) return sim_clean();

                        flag = PyTuple_SetItem(l1, 0, l2); /* Steals reference to l2 */
                        l2 = NULL; flt = NULL;

                        l2 = PyTuple_New(model->n_sens);
                        if (l2 == NULL) return sim_clean();

                        flt = PyFloat_FromDouble(model->S0_V_INa);
                        if (flt == NULL) return sim_clean();
                        flag = PyTuple_SetItem(l2, 0, flt); /* Steals reference to flt */
                        if (flag < 0) return sim_clean();
                        flt = PyFloat_FromDouble(model->S1_V_INa);
                        if (flt == NULL) return sim_clean();
                        flag = PyTuple_SetItem(l2, 1, flt); /* Steals reference to flt */
                        if (flag < 0) return sim_clean();
                        flt = PyFloat_FromDouble(model->S2_V_INa);
                        if (flt == NULL) return sim_clean();
                        flag = PyTuple_SetItem(l2, 2, flt); /* Steals reference to flt */
                        if (flag < 0) return sim_clean();

                        flag = PyTuple_SetItem(l1, 1, l2); /* Steals reference to l2 */
                        l2 = NULL; flt = NULL;

                        flag = PyList_Append(sens_list, l1);
                        Py_XDECREF(l1); l1 = NULL;
                        if (flag < 0) return sim_clean();

                    }

                    /* Get next logging point */
                    if (log_interval > 0) {
                        /* Periodic logging */
                        ilog++;
                        tlog = tmin + (double)ilog * log_interval;
                        if (ilog == 0) {
                            /* Unsigned int wraps around instead of overflowing, becomes zero again */
                            PyErr_SetString(PyExc_Exception, "Overflow in logged step count: Simulation too long!");
                            return sim_clean();
                        }
                    } else {
                        /* Point-list logging */
                        /* Read next log point off the list */
                        if (ilog < PyList_Size(log_times)) {
                            flt = PyList_GetItem(log_times, ilog); /* Borrowed */
                            if (!PyFloat_Check(flt)) {
                                PyErr_SetString(PyExc_Exception, "Entries in 'log_times' must be floats.");
                                return sim_clean();
                            }
                            tlog = PyFloat_AsDouble(flt);
                            ilog++;
                            flt = NULL;
                        } else {
                            tlog = tmax + 1;
                        }
                    }
                }
            }

            /* Event-based pacing */

            /* At this point we have logged everything _before_ engine_time, so
               it's safe to update the pacing mechanism. */
            if (epacing != NULL) {
                flag_epacing = ESys_AdvanceTime(epacing, engine_time);
                if (flag_epacing != ESys_OK) { ESys_SetPyErr(flag_epacing); return sim_clean(); }
                tnext = ESys_GetNextTime(epacing, NULL);
                engine_pace = ESys_GetLevel(epacing, NULL);
                tnext = (tnext < tmax) ? tnext : tmax;
            }

            /* Dynamic logging: Log every visited point */
            if (dynamic_logging) {

                /* Benchmarking? Then set engine_realtime */
                if (benchtime != Py_None) {
                    flt = PyObject_CallFunction(benchtime, "");
                    if (!PyFloat_Check(flt)) {
                        Py_XDECREF(flt); flt = NULL;
                        PyErr_SetString(PyExc_Exception, "Call to benchmark time function didn't return float.");
                        return sim_clean();
                    }
                    engine_realtime = PyFloat_AsDouble(flt) - engine_starttime;
                    Py_DECREF(flt); flt = NULL;
                }

                /* Ensure the logged values are correct for the new time t */
                if (log_deriv || log_inter) {
                    /* If logging derivatives or intermediaries, calculate the
                       values for the current time. */
                    /*TODO:REPLACE */
                    rhs(engine_time, y, dy_log, udata);
                } else if (log_bound) {
                    /* Logging bounds but not derivs or inters: No need to run
                       full rhs, just update bound variables */
                    Model_UpdateBoundVariables(
                        model, engine_time, engine_pace,
                        engine_realtime, engine_evaluations);
                }

                /* Write to log */
                for (i=0; i<n_vars; i++) {
                    flt = PyFloat_FromDouble(*vars[i]);
                    ret = PyObject_CallMethodObjArgs(logs[i], list_update_str, flt, NULL);
                    Py_DECREF(flt); flt = NULL;
                    Py_XDECREF(ret);
                    if (ret == NULL) {
                        PyErr_SetString(PyExc_Exception, "Call to append() failed on logging list.");
                        return sim_clean();
                    }
                    ret = NULL;
                }

                if (model->n_sens) {

                    /* Calculate sensitivities to output */
                    calculate_sensitivity_outputs(engine_time, y, dy_log, sy, sdy_log, udata);


                    /* Write sensitivity matrix to log */
                    l1 = PyTuple_New(2);
                    if (l1 == NULL) return sim_clean();

                    l2 = PyTuple_New(model->n_sens);
                    if (l2 == NULL) return sim_clean();

                    flt = PyFloat_FromDouble(model->S0_Y_V);
                    if (flt == NULL) return sim_clean();
                    flag = PyTuple_SetItem(l2, 0, flt); /* Steals reference to flt */
                    if (flag < 0) return sim_clean();
                    flt = PyFloat_FromDouble(model->S1_Y_V);
                    if (flt == NULL) return sim_clean();
                    flag = PyTuple_SetItem(l2, 1, flt); /* Steals reference to flt */
                    if (flag < 0) return sim_clean();
                    flt = PyFloat_FromDouble(model->S2_Y_V);
                    if (flt == NULL) return sim_clean();
                    flag = PyTuple_SetItem(l2, 2, flt); /* Steals reference to flt */
                    if (flag < 0) return sim_clean();

                    flag = PyTuple_SetItem(l1, 0, l2); /* Steals reference to l2 */
                    l2 = NULL; flt = NULL;

                    l2 = PyTuple_New(model->n_sens);
                    if (l2 == NULL) return sim_clean();

                    flt = PyFloat_FromDouble(model->S0_V_INa);
                    if (flt == NULL) return sim_clean();
                    flag = PyTuple_SetItem(l2, 0, flt); /* Steals reference to flt */
                    if (flag < 0) return sim_clean();
                    flt = PyFloat_FromDouble(model->S1_V_INa);
                    if (flt == NULL) return sim_clean();
                    flag = PyTuple_SetItem(l2, 1, flt); /* Steals reference to flt */
                    if (flag < 0) return sim_clean();
                    flt = PyFloat_FromDouble(model->S2_V_INa);
                    if (flt == NULL) return sim_clean();
                    flag = PyTuple_SetItem(l2, 2, flt); /* Steals reference to flt */
                    if (flag < 0) return sim_clean();

                    flag = PyTuple_SetItem(l1, 1, l2); /* Steals reference to l2 */
                    l2 = NULL; flt = NULL;

                    flag = PyList_Append(sens_list, l1);
                    Py_XDECREF(l1); l1 = NULL;
                    if (flag < 0) return sim_clean();

                }

            }

            /* Reinitialize CVODE if needed (cvode-mode only) */
            if (model->is_ode && flag_reinit) {
                flag_reinit = 0;
                /* Re-init */
                flag_cvode = CVodeReInit(cvode_mem, engine_time, y);
                if (check_cvode_flag(&flag_cvode, "CVodeReInit", 1)) return sim_clean();
                flag_cvode = CVodeSensReInit(cvode_mem, CV_SIMULTANEOUS, sy);
                if (check_cvode_flag(&flag_cvode, "CVodeSensReInit", 1)) return sim_clean();
            }
        }

        /* Check if we're finished */
        if (ESys_eq(engine_time, tmax)) engine_time = tmax;
        if (engine_time >= tmax) break;

        /* Perform any Python signal handling */
        if (PyErr_CheckSignals() != 0) {
            /* Exception (e.g. timeout or keyboard interrupt) occurred?
               Then cancel everything! */
            return sim_clean();
        }

        /* Report back to python after every x steps */
        steps_taken++;
        if (steps_taken >= 100) {
            return PyFloat_FromDouble(engine_time);
        }
    }

    /* Set final state */
    for (i=0; i<model->n_states; i++) {
        PyList_SetItem(state_out, i, PyFloat_FromDouble(NV_Ith_S(y, i)));
        /* PyList_SetItem steals a reference: no need to decref the double! */
    }

    /* Set bound variable values */
    PyList_SetItem(inputs, 0, PyFloat_FromDouble(engine_time));
    PyList_SetItem(inputs, 1, PyFloat_FromDouble(engine_pace));
    PyList_SetItem(inputs, 2, PyFloat_FromDouble(engine_realtime));
    PyList_SetItem(inputs, 3, PyFloat_FromDouble(engine_evaluations));

    sim_clean();    /* Ignore return value */
    return PyFloat_FromDouble(engine_time);
}

/*
 * Evaluates the state derivatives at the given state
 */
static PyObject*
sim_eval_derivatives(PyObject *self, PyObject *args)
{
    /* Declare variables here for C89 compatibility */
    int i;
    int success;
    int iState;
    double time_in;
    double pace_in;
    char errstr[200];
    Model_Flag flag_model;
    PyObject *state;
    PyObject *deriv;
    PyObject *flt;
    Model m;
    N_Vector y;
    N_Vector dy;
    /* TODO: CREATE OWN UDATA HERE */

    /* Start */
    success = 0;

    /* Check input arguments */
    if (!PyArg_ParseTuple(args, "OOdd", &state, &deriv, &time_in, &pace_in)) {
        PyErr_SetString(PyExc_Exception, "Expecting sequence arguments 'y' and 'dy' followed by floats 'time' and 'pace'.");
        /* Nothing allocated yet, no pyobjects _created_, return directly */
        return 0;
    }
    if (!PySequence_Check(state)) {
        PyErr_SetString(PyExc_Exception, "First argument must support the sequence interface.");
        return 0;
    }
    if (!PySequence_Check(deriv)) {
        PyErr_SetString(PyExc_Exception, "Second argument must support the sequence interface.");
        return 0;
    }

    /* From this point on, no more direct returning: use goto error */
    m = NULL;
    y = NULL;      /* A cvode SERIAL vector */
    dy = NULL;     /* A cvode SERIAL vector */
    udata = NULL;  /* The user data, containing parameter values */
    #if MYOKIT_SUNDIALS_VERSION >= 30000
    sundense_matrix = NULL;     /* A matrix for linear solving */
    sundense_solver = NULL;     /* A linear solver */
    #endif

    /* Temporary object: decref before re-using for another var :) */
    /* (Unless you get them using PyList_GetItem...) */
    flt = NULL;   /* PyFloat */

    /* Create model */
    m = Model_Create(&flag_model);
    if (flag_model != Model_OK) {
        Model_SetPyErr(flag_model);
        goto error;
    }

    /* Create state vectors */
    y = N_VNew_Serial(m->n_states);
    if (check_cvode_flag((void*)y, "N_VNew_Serial", 0)) {
        PyErr_SetString(PyExc_Exception, "Failed to create state vector.");
        goto error;
    }
    dy = N_VNew_Serial(m->n_states);
    if (check_cvode_flag((void*)dy, "N_VNew_Serial", 0)) {
        PyErr_SetString(PyExc_Exception, "Failed to create state derivatives vector.");
        goto error;
    }

    /* Set calculated constants (not including sensitivity parameters) */
    Model_UpdateLiteralDerivedVariables(model);

    /* Set initial values */
    for (iState = 0; iState < m->n_states; iState++) {
        flt = PySequence_GetItem(state, iState); /* Remember to decref! */
        if (!PyFloat_Check(flt)) {
            Py_XDECREF(flt); flt = NULL;
            sprintf(errstr, "Item %d in state vector is not a float.", iState);
            PyErr_SetString(PyExc_Exception, errstr);
            goto error;
        }
        NV_Ith_S(y, iState) = PyFloat_AsDouble(flt);
        Py_DECREF(flt);
    }
    flt = NULL;

    /* Create and fill vector with parameter/initial condition values */
    udata = (UserData)malloc(sizeof *udata);
    if (udata == 0) {
        PyErr_SetString(PyExc_Exception, "Unable to allocate space to store parameter values.");
        return sim_clean();
    }
    udata->p[0] = model->P_gNa;
    udata->p[1] = model->P_gCa;
    udata->p[2] = NV_Ith_S(y, 7);


    /* Set simulation time and pacing variable */
    engine_time = time_in;
    engine_pace = pace_in;

    /* Evaluate derivatives */
    rhs(engine_time, y, dy, udata);

    /* Set output values */
    for (i=0; i<model->n_states; i++) {
        flt = PyFloat_FromDouble(NV_Ith_S(dy, i));
        if (flt == NULL) {
            PyErr_SetString(PyExc_Exception, "Unable to create float.");
            goto error;
        }
        PySequence_SetItem(deriv, i, flt);
        Py_DECREF(flt);
    }
    flt = NULL;

    /* Finished succesfully, free memory and return */
    success = 1;
error:
    /* Free CVODE space */
    if (y != NULL) { N_VDestroy_Serial(y); }
    if (dy != NULL) { N_VDestroy_Serial(dy); }

    /* Free udata space */
    free(udata->p);
    free(udata);

    /* Free model space */
    Model_Destroy(m);

    /* Return */
    if (success) {
        Py_RETURN_NONE;
    } else {
        return 0;
    }
}

/*
 * Returns the number of steps taken in the last simulation
 */
static PyObject*
sim_steps(PyObject *self, PyObject *args)
{
    return PyLong_FromLong(engine_steps);
}

/*
 * Returns the number of rhs evaluations performed during the last simulation
 */
static PyObject*
sim_evals(PyObject *self, PyObject *args)
{
    return PyLong_FromLong(engine_evaluations);
}

/*
 * Methods in this module
 */
static PyMethodDef SimMethods[] = {
    {"sim_init", sim_init, METH_VARARGS, "Initialize the simulation."},
    {"sim_step", sim_step, METH_VARARGS, "Perform the next step in the simulation."},
    {"sim_clean", py_sim_clean, METH_VARARGS, "Clean up after an aborted simulation."},
    {"eval_derivatives", sim_eval_derivatives, METH_VARARGS, "Evaluate the state derivatives."},
    {"set_tolerance", sim_set_tolerance, METH_VARARGS, "Set the absolute and relative solver tolerance."},
    {"set_max_step_size", sim_set_max_step_size, METH_VARARGS, "Set the maximum solver step size (0 for none)."},
    {"set_min_step_size", sim_set_min_step_size, METH_VARARGS, "Set the minimum solver step size (0 for none)."},
    {"number_of_steps", sim_steps, METH_VARARGS, "Returns the number of steps taken in the last simulation."},
    {"number_of_evaluations", sim_evals, METH_VARARGS, "Returns the number of rhs evaluations performed during the last simulation."},
    {NULL},
};

/*
 * Module definition
 */
#if PY_MAJOR_VERSION >= 3

    static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "myokit_sim_1_4384500282905115256",       /* m_name */
        "Generated CVODESim module",/* m_doc */
        -1,                         /* m_size */
        SimMethods,                 /* m_methods */
        NULL,                       /* m_reload */
        NULL,                       /* m_traverse */
        NULL,                       /* m_clear */
        NULL,                       /* m_free */
    };

    PyMODINIT_FUNC PyInit_myokit_sim_1_4384500282905115256(void) {
        return PyModule_Create(&moduledef);
    }

#else

    PyMODINIT_FUNC
    initmyokit_sim_1_4384500282905115256(void) {
        (void) Py_InitModule("myokit_sim_1_4384500282905115256", SimMethods);
    }

#endif

