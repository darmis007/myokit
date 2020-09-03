
/*
 * Right-hand-side function of the model's sensitivity equations.
 *
 *  int Ns          The number of sensitivities
 *  realtype t      Current time
 *  N_Vector y      The current state values
 *  N_Vector ydot   The current state value derivatives
 *  N_Vector* yS    The current values of the sensitivity vectors (see below)
 *  N_Vector* ySdot Space to store the calculated sensitivities of the derivatives
 *  void* user_data Extra data (contains the sensitivity parameter values)
 *  N_Vector tmp1   A length-N N_Vector for use as temporary storage
 *  N_Vector tmp2   A length-N N_Vector for use as temporary storage
 *
 * yS[i] is a vector containing s_i = dy/dp_i
 *  Each entry in yS is for one parameter/initial value
 *  And each entry in yS[i] is for a state
 *
 * ySdot[i] is a vector containing (df/dy) * s_i + df/dp_i
 *  Each entry in ySdot is for one parameter/initial value
 *  And each entry in ySdot[i] is for a state
 */
static int
rhs1(int Ns, realtype t, N_Vector y, N_Vector ydot, N_Vector *yS, N_Vector *ySdot,
   void *user_data, N_Vector tmp1, N_Vector tmp2)
{
    FSys_Flag flag_fpacing;
    UserData fdata;

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

<?
'''
for label, eqs in equations.items():
    need_label = True
    for eq in eqs.equations():
        var = eq.lhs.var()

        # Skip parameters
        if var in parameters:
            continue

        # Skip constants that don't depend on parameters
        if var.is_constant() and len(eq.rhs.references().intersection(parameter_name_set)) == 0:
            continue

        # Print label for component
        if need_label:
            print(tab + '/* ' + label + ' */')
            need_label = False

        # Print equation
        try:
            print(tab + v(var) + ' = ' + bound_variables[var] + ';')
        except KeyError:
            print(tab + w.eq(eq) + ';')

    if not need_label:
        print(tab)

for i, eqs in enumerate(sensitivity_equations):
    print('')
    print(tab + '/* Sensitivities for p[' + str(i) + '] */')
    for eq in eqs:
        print(tab + w.eq(eq) + ';')

print('')
print(tab + '/* Calcuate (df/dy)s_i(t) + df/dp_i */')
for i, expr in enumerate(indep_expr):
    #
    # ysDot[i] = df/dpi + df/dy * si
    #
    # Where df/dy is a matrix of size (n_states, n_states)
    #
    # For each sensitivity i, for each row j, we get
    #
    # ysDot[i][j] = df[j]/dp[i] + Sum_k df[j]/dy[k] * s[i][k]
    #             = df[j]/dp[i] + Sum_k df[j]/dy[k] * dy[k]/dp[i]
    #
    for j, yj in enumerate(model.states()):
        # Initialise as: yDot[i][j] = df[j]/dp[i]
        #                           = d/dp[i] dot(y[j])
        dfj_dpi = myokit.PartialDerivative(myokit.Derivative(myokit.Name(yj)), expr)
        print(tab + 'NV_Ith_S(ySdot[' + str(i) + '], ' + str(j) + ') = ' + v(dfj_dpi), end=None)
        for k, yk in enumerate(model.states()):
            dfj_dyk = v(myokit.PartialDerivative(myokit.Derivative(myokit.Name(yj)), myokit.Name(yk)))
            dyk_dpi = 'NV_Ith_S(yS[' + str(i) + '], ' + str(k) + ')'
            print(' + ' + dfj_dyk + ' * ' + dyk_dpi, end=None)
        print(';')
'''
?>
    return 0;
}

