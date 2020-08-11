#
# Converts SBML to Myokit expressions, using an ElementTree implementation.
#
# Only partial SBML support (based on SBML level 3 version 2) is provided.
# The SBML file format specifications can be found here
# http://sbml.org/Documents/Specifications.
#
# This file is part of Myokit.
# See http://myokit.org for copyright, sharing, and licensing details.
#
from __future__ import absolute_import, division
from __future__ import print_function, unicode_literals

import collections
import re
import warnings

import xml.etree.ElementTree as ET

import myokit
import myokit.units
import myokit.formats
import myokit.formats.mathml
import myokit.formats.html
import myokit.formats.xml


# Namespaces
NS_MATHML = 'http://www.w3.org/1998/Math/MathML'
NS_SBML_3_2 = 'http://www.sbml.org/sbml/level3/version2/core'


# Regex for id checking
_re_id = re.compile(r'^[a-zA-Z_]+[a-zA-Z0-9_]*$')


def _parse_mathml(expr, name_generator):
    """
    Parses the MathML expression ``expr`` and wraps any errors in an SBMLError.
    """
    def name(x, y):
        try:
            return name_generator(x, y)
        except KeyError:
            raise SBMLError(
                'Unknown or inaccessible symbol in MathML: "' + str(x) + '".')

    try:
        return myokit.formats.mathml.parse_mathml_etree(
            expr, name, lambda x, y: myokit.Number(x))
    except myokit.formats.mathml.MathMLError as e:
        raise SBMLError('Unable to parse MathML: ' + str(e))

    # Note: In any place we parse mathematics, we could also check that
    # there's nothing after the mathematics. E.g.
    #   <math xmlns="..">
    #    <cn>1</cn>
    #    <cn>2</cn>
    #   </math>
    # is currently allowed because only the first statement is read.
    # But for now we're being lenient.


class SBMLError(Exception):
    """Raised if something goes wrong when working with an SBML model."""


class SBMLParsingError(myokit.ImportError):
    """
    Thrown if an error occurs when parsing SBML.

    The argument ``element`` can be used to pass in an element that caused the
    error.
    """
    def __init__(self, message, element=None):
        if element is not None:
            try:
                line = str(element.sourceline)
                message = 'Error on line ' + line + '. ' + message
            except AttributeError:  # pragma: no cover
                pass
        super(SBMLParsingError, self).__init__(message)


class SBMLParser(object):
    """
    Parses SBML models, creating an SBML Model structure that can be converted
    to a :class:`myokit.Model` object.

    General notes:

    - Component and variable names are created from ``id`` attributes, and
      ``name`` attributes are ignored. (Names are not typically required,
      while ids are; ids are almost myokit identifiers already (except that
      they can start with an underscore); and ids are guaranteed to be
      unique).

    Support notes:

    - SBML older than Level 2 version 1 is unlikely to work.
    - Algebraic rules are not supported.
    - Constraints are not supported.
    - Events are not supported.
    - Function definitions are not supported.
    - Local parameters in reactions are not supported.
    - Units "celsius" are not supported.
    - Since SBML does not support units in equations, variable units will be
      set, but RHS units will not.

    """

    def _check_required_attributes(self, element, *attributes):
        """
        Checks that an ``elementtree`` element has all the required attributes.
        """
        for attribute in attributes:
            if element.get(attribute) is None:
                raise SBMLParsingError(
                    'Element ' + self._tag(element) + ' is missing required'
                    ' attribute "' + str(attribute) + '".', element)

    def parse_file(self, path):
        """
        Parses the SBML file at ``path`` and returns a :class:`myokit.Model`.
        """
        # Read file
        try:
            tree = ET.parse(path)
        except Exception as e:
            raise SBMLParsingError('Unable to parse XML: ' + str(e))

        # Parse content
        return self.parse(tree.getroot())

    def parse_string(self, text):
        """
        Parses the SBML XML in the string ``text`` and returns a
        :class:`myokit.Model`.
        """
        # Read string
        try:
            root = ET.fromstring(text)
        except Exception as e:
            raise SBMLParsingError('Unable to parse XML: ' + str(e))

        # Parse content
        return self.parse(root)

    def parse(self, element):
        """
        Parses an SBML document rooted in the given ``ElementTree`` element and
        returns a :class:`myokit.Model`.
        """
        # Supported namespaces
        # Other namespaces are allowed, but might not work.
        supported = (
            NS_SBML_3_2,
        )

        # Check whether document declares a supported namespace
        self._ns = myokit.formats.xml.split(element.tag)[0]
        try:
            if self._ns not in supported:
                warnings.warn(
                    'Unknown SBML namespace ' + str(self._ns) + '. This'
                    ' version of SBML may not be supported.')

            # Get model and parse
            model = element.find(self._path('model'))
            if model is None:
                raise SBMLParsingError('Model element not found.', element)
            return self._parse_model(model)

        finally:
            # Remove all references to temporary state
            del(self._ns)

    def _parse_compartment(self, element, model):
        """Parses a ``compartment`` element."""

        # Check required attributes
        self._check_required_attributes(element, 'id')
        # The attribute 'constant' is also required in level 3, but optional in
        # level 2 (default true). If constant is false then the compartment's
        # size can be changed.

        # Note: Any place we check required attributes, we could also check for
        # allowed attributes (e.g. errors for unrecognised ones). Bit tricky
        # for multiple SBML versions. For now we're being lenient.

        sid = element.get('id')
        try:
            # Create compartment
            compartment = model.add_compartment(sid)

            # Get units, if given
            units = element.get('units')
            if units is not None:
                try:
                    units = model.unit(units)
                except KeyError:
                    raise SBMLParsingError(
                        'Unknown units "' + units + '".', element)
                compartment.set_size_units(units)

            # Get spatial dimensions, if given
            dim = element.get('spatialDimensions')
            if dim is not None:
                try:
                    dim = float(dim)
                except ValueError:
                    raise SBMLParsingError(
                        'Unable to convert spatial dimensions value "'
                        + str(dim) + '" to float.', element)
                compartment.set_spatial_dimensions(dim)

            # Get compartment size, if set
            size = element.get('size')
            if size is not None:
                try:
                    size = float(size)
                except ValueError:
                    raise SBMLParsingError(
                        'Unable to convert size "' + str(size) + '" to float.',
                        element)

                # Ignore units in RHS
                compartment.set_initial_value(myokit.Number(size))

        except SBMLError as e:
            raise SBMLParsingError(
                'Unable to parse compartment "' + sid + '": ' + str(e),
                element)

    def _parse_initial_assignment(self, element, model):
        """Parses an initial assignment."""

        # Check required attributes
        self._check_required_attributes(element, 'symbol')

        # Check that symbol refers to an assignable
        var = element.get('symbol')
        try:
            var = model.assignable(var)
        except KeyError:
            raise SBMLParsingError(
                'Unable to parse rule: The given SId "' + str(var) + '" does'
                ' not refer to a Compartment, Species, SpeciesReference, or'
                ' Parameter in this model.')

        try:
            # Find mathematics
            expr = element.find('{' + NS_MATHML + '}math')
            if expr is None:
                # This is allowed, so that SBML extensions can do mysterious
                # things.
                warnings.warn(
                    'Initial assignment does not define any mathematics')
                return

            # Parse rule and apply
            var.set_initial_value(_parse_mathml(
                expr,
                lambda x, y: myokit.Name(model.assignable(x)),
            ))

            # Indicate for species that the units of the inital value are
            # correct (initialAssignments are meant to assign the value
            # in the correct units amount/concentration)
            if isinstance(var, Species):
                var.set_units_initial_value(True)

        except SBMLError as e:
            raise SBMLParsingError(
                'Unable to parse initial assignment: ' + str(e), element)

    def _parse_kinetic_law(self, element, model, reaction):
        """
        Parses the ``kinetic_law`` element used in reactions and returns the
        resulting :class:`myokit.Expression``.

        Returns ``None`` if no expression is found.
        """

        # Local parameters are not supported
        x = element.find(self._path('listOfLocalParameters', 'localParameter'))
        if x is not None:
            raise SBMLParsingError('Local parameters are not supported.', x)

        # Kinetic law is not allowed use species unless they're defined in the
        # reaction.
        def name(x, y):
            obj = model.assignable(x)
            if isinstance(obj, Species):
                obj = reaction.species(x)
            return myokit.Name(obj)

        # Find and parse expression
        expr = element.find('{' + NS_MATHML + '}math')
        if expr is None:
            return None
        try:
            return _parse_mathml(expr, name)
        except SBMLError as e:
            raise SBMLParsingError(
                'Unable to parse kinetic law: ' + str(e), element)

    def _parse_model(self, element):
        """Parses a ``model`` element in an SBML document."""

        # Retrieve or create a model name.
        # SBML Models can have an optional name attribute (user friendly name)
        # or an optional id attribute (not necessarily user friendly) or can
        # have no name at all.
        name = element.get('name')
        if name is None:
            name = element.get('id')
        if name is None:
            name = 'Imported SBML model'

        # Create SBML model
        model = Model(name)

        # Algebraic rules are not supported
        x = element.find(self._path('listOfRules', 'algebraicRule'))
        if x is not None:
            raise SBMLParsingError(
                'Algebraic rules are not supported.', x)

        # Constraints are not supported, but will be ignored
        x = element.find(self._path('listOfConstraints', 'constraint'))
        if x is not None:
            warnings.warn('Ignoring SBML constraints.')

        # Events are not supported, but will be ignored
        x = element.find(self._path('listOfEvents', 'event'))
        if x is not None:
            warnings.warn('Ignoring SBML events.')

        # Function definitions are not supported.
        x = element.find(self._path(
            'listOfFunctionDefinitions', 'functionDefinition'))
        if x is not None:
            raise SBMLParsingError(
                'Function definitions are not supported.', x)

        # Parse notes element (converts HTML to plain text)
        notes = element.find(self._path('notes'))
        if notes is not None:
            self._parse_notes(notes, model)

        # Parse unit definitions
        for unit in element.findall(self._path(
                'listOfUnitDefinitions', 'unitDefinition')):
            self._parse_unit_definition(unit, model)

        # Parse model unit attributes
        try:
            # Default units for compartment sizes
            units = element.get('lengthUnits')
            if units is not None:
                model.set_length_units(model.unit(units))
            units = element.get('areaUnits')
            if units is not None:
                model.set_area_units(model.unit(units))
            units = element.get('volumeUnits')
            if units is not None:
                model.set_volume_units(model.unit(units))

            # Default units for species amounts (not concentrations)
            units = element.get('substanceUnits')
            if units is not None:
                model.set_substance_units(model.unit(units))

            # Default units for reactions and rates.
            # The KineticLaw equations are in units of extentUnits / timeUnits.
            units = element.get('extentUnits')
            if units is not None:
                model.set_extent_units(model.unit(units))
            units = element.get('timeUnits')
            if units is not None:
                model.set_time_units(model.unit(units))

        except KeyError as e:
            raise SBMLParsingError(
                'Unknown units "' + str(e.args[0]) + '".', element)
        except SBMLError as e:
            raise SBMLParsingError(
                'Error parsing model element: ' + str(e), element)

        # Parse parameters
        for parameter in element.findall(self._path(
                'listOfParameters', 'parameter')):
            self._parse_parameter(parameter, model)

        # Set global conversion factor, if provided (must refer to a parameter)
        factor = element.get('conversionFactor')
        if factor is not None:
            try:
                factor = model.parameter(factor)
            except KeyError:
                raise SBMLParsingError(
                    'Model conversion factor "' + str(factor) + '" does not'
                    ' refer to a parameter SId.', element)
            model.set_conversion_factor(factor)

        # Parse compartments
        for compartment in element.findall(self._path(
                'listOfCompartments', 'compartment')):
            self._parse_compartment(compartment, model)

        # Parse species
        for species in element.findall(self._path('listOfSpecies', 'species')):
            self._parse_species(species, model)

        # Parse reactions
        for reaction in element.findall(self._path(
                'listOfReactions', 'reaction')):
            self._parse_reaction(reaction, model)

        # Parse initial assignments
        for assignment in element.findall(self._path(
                'listOfInitialAssignments', 'initialAssignment')):
            self._parse_initial_assignment(assignment, model)

        # Parse assignment rules
        for rule in element.findall(self._path(
                'listOfRules', 'assignmentRule')):
            self._parse_rule(rule, model, False)

        # Parse rate rules
        for rule in element.findall(self._path('listOfRules', 'rateRule')):
            self._parse_rule(rule, model, True)

        return model

    def _parse_notes(self, element, model):
        """Parses a model's ``notes`` element, converting it to plain text."""

        notes = ET.tostring(element).decode()
        notes = myokit.formats.html.html2ascii(notes, width=75)
        if notes:
            model.set_notes(notes)

    def _parse_parameter(self, element, model):
        """Parses a ``parameter`` element."""

        # Check required attributes
        self._check_required_attributes(element, 'id')
        # The attribute 'constant' is also required in level 3, but optional in
        # level 2 (default false)

        # Create
        try:
            parameter = model.add_parameter(element.get('id'))

            # Get optional units
            unit = element.get('units')
            if unit is not None:
                try:
                    unit = model.unit(unit)
                except KeyError:
                    raise SBMLParsingError(
                        'Unknown units "' + unit + '".', element)
                parameter.set_units(unit)

            # Get optional initial value
            value = element.get('value')
            if value is not None:
                try:
                    value = float(value)
                except ValueError:
                    raise SBMLParsingError(
                        'Unable to convert parameter value "' + str(value)
                        + '" to float.', element)

                # Ignore unit
                parameter.set_initial_value(myokit.Number(value))

        except SBMLError as e:
            raise SBMLParsingError(
                'Unable to parse ' + self._tag(element) + ': ' + str(e),
                element)

    def _parse_reaction(self, element, model):
        """Parses a ``reaction`` element."""

        # Check required attributes
        self._check_required_attributes(element, 'id')
        # The attribute 'reversible' is also required, but doesn't affect the
        # maths, so we ignore it here.
        # The attribute 'compartment' is not required, and doesn't affect the
        # maths, so we ignore it here.

        # The 'fast' attribute was removed in L3V2 and is not supported here
        if element.get('fast') == 'true':
            raise SBMLParsingError(
                'Reactions with fast="true" are not supported.', element)

        # Create
        try:
            reaction = model.add_reaction(element.get('id'))
        except SBMLError as e:
            raise SBMLParsingError(
                'Unable to parse ' + self._tag(element) + ': ' + str(e),
                element)

        # Add reactants, products, modifiers
        have_reactant_or_product = False
        for x in element.findall(self._path(
                'listOfReactants', 'speciesReference')):
            self._parse_species_reference(
                x, model, reaction, SpeciesReference.REACTANT)
            have_reactant_or_product = True

        for x in element.findall(self._path(
                'listOfProducts', 'speciesReference')):
            self._parse_species_reference(
                x, model, reaction, SpeciesReference.PRODUCT)
            have_reactant_or_product = True

        for x in element.findall(self._path(
                'listOfModifiers', 'modifierSpeciesReference')):
            self._parse_species_reference(
                x, model, reaction, SpeciesReference.MODIFIER)

        # Raise error if neither reactants nor products is populated
        if not have_reactant_or_product:
            raise SBMLParsingError(
                'Reaction must have at least one reactant or product.',
                element)

        # Parse kinetic law
        kinetic_law = element.find(self._path('kineticLaw'))
        if kinetic_law is not None:
            kinetic_law = self._parse_kinetic_law(kinetic_law, model, reaction)
        if kinetic_law is None:
            warnings.warn(
                'No kinetic law set for reaction "' + reaction.sid() + '".')
        else:
            reaction.set_kinetic_law(kinetic_law)

    def _parse_rule(self, element, model, rate=False):
        """Parses an assignment or rate rule."""

        # Check required attributes
        self._check_required_attributes(element, 'variable')

        # Check that variable refers to an assignable
        var = element.get('variable')
        try:
            var = model.assignable(var)
        except KeyError:
            raise SBMLParsingError(
                'Unable to parse rule: The given SId "' + str(var) + '" does'
                ' not refer to a Compartment, Species, SpeciesReference, or'
                ' Parameter in this model.')

        try:
            # Check that this assignable can be changed with assignment rules
            if isinstance(var, Species):
                if not var.boundary():
                    raise SBMLParsingError(
                        'Assignment or rate rule set for species that is'
                        ' created or destroyed in a reaction.', element)
            # TODO: Could also gather 'constant' attributes for all, and check
            # here. But maybe we're not all that interested in validation...

            # Find mathematics
            expr = element.find('{' + NS_MATHML + '}math')
            if expr is None:
                # This is allowed, so that SBML extensions can do mysterious
                # things.
                warnings.warn('Rule does not define any mathematics')
                return

            # Parse rule and apply
            var.set_value(_parse_mathml(
                expr,
                lambda x, y: myokit.Name(model.assignable(x)),
            ), rate)

        except SBMLError as e:
            raise SBMLParsingError('Unable to parse rule: ' + str(e), element)

    def _parse_species(self, element, model):
        """
        Adds species to references compartment in model.
        """
        # Check required attributes (note: constant is also required in L3)
        self._check_required_attributes(element, 'id', 'compartment')
        # The attributes 'hasOnlySubstaneUnits', 'boundaryCondition', and
        # 'constant' are required in level 3, but optional in level 2 (default
        # is false for all three).

        # Check if it's an amount or a concentration
        amount = element.get('hasOnlySubstanceUnits', 'false') == 'true'

        # Check if constant, and if at a reaction boundary
        constant = element.get('constant', 'false') == 'true'
        boundary = element.get('boundaryCondition', 'false') == 'true'

        # Note: In lines like the above we could raise an error if the value
        # isn't 'true' or false', but for now we're being lenient.

        # Get compartment
        compartment = element.get('compartment')
        try:
            compartment = model.compartment(compartment)
        except KeyError:
            raise SBMLParsingError(
                'Unknown compartment "' + compartment + '"', element)

        # Create
        try:
            species = model.add_species(
                compartment, element.get('id'), amount, constant, boundary)

            # Set units, if provided
            units = element.get('substanceUnits')
            if units is not None:
                try:
                    units = model.unit(units)
                except KeyError:
                    raise SBMLParsingError(
                        'Unknown units "' + units + '"', element)

                # Set substance units, not concentration units
                species.set_substance_units(units)

            # Check that not both initialAmount and initialConcentration are
            # provided
            if (element.get('initialAmount') is not None) and (
                    element.get('initialConcentration') is not None):
                raise SBMLParsingError(
                    'Species cannot set both an initialAmount and an'
                    ' initialConcentration.')

            # Get initial amount
            value = element.get('initialAmount')

            # Indicate whether units of initial value are correct
            species.set_units_initial_value(species.amount())

            # If initial amount is not provided, get initial concentration
            if value is None:
                value = element.get('initialConcentration')

                # Indicate whether units of initial value are correct
                species.set_units_initial_value(not species.amount())

            # Set initial value
            if value is not None:
                try:
                    value = float(value)
                except ValueError:
                    raise SBMLParsingError(
                        'Unable to convert initial species value to float "'
                        + str(value) + '".', element)

                # Ignore units in equations
                species.set_initial_value(myokit.Number(value))

            # Set conversion factor if provided (must refer to a parameter)
            factor = element.get('conversionFactor')
            if factor is not None:
                try:
                    factor = model.parameter(factor)
                except KeyError:
                    raise SBMLParsingError(
                        'Unknown parameter "' + str(factor) + '" set as'
                        ' conversion factor.', element)
                species.set_conversion_factor(factor)

        except SBMLError as e:
            raise SBMLParsingError(
                'Unable to parse ' + self._tag(element) + ': ' + str(e),
                element)

    def _parse_species_reference(self, element, model, reaction, ref_type):
        """Parses a ``speciesReference`` element inside a reaction."""

        # Check required attributes
        # Note: Level 3 also has a required attribute 'constant'
        self._check_required_attributes(element, 'species')

        # Find species
        species = element.get('species')
        try:
            species = model.species(species)
        except KeyError:
            raise SBMLParsingError(
                'Reference to unknown species "' + species + '"', element)

        # Get optional SId
        sid = element.get('id')

        try:
            # Create
            if ref_type == SpeciesReference.REACTANT:
                reference = reaction.add_reactant(species, sid)
            elif ref_type == SpeciesReference.PRODUCT:
                reference = reaction.add_product(species, sid)
            else:
                # Modifier references require no further parsing
                return reaction.add_modifier(species, sid)

            # Get optional stoichiometry
            stoichiometry = element.get('stoichiometry')
            if stoichiometry is not None:
                try:
                    stoichiometry = float(stoichiometry)
                except ValueError:
                    raise SBMLParsingError(
                        'Unable to convert stoichiometry value "'
                        + str(stoichiometry) + '" to float.', element)
                reference.set_initial_value(myokit.Number(stoichiometry))

        except SBMLError as e:
            raise SBMLParsingError(
                'Unable to parse species reference: ' + str(e), element)

        return reference

    def _parse_unit_definition(self, element, model):
        """Parses a ``unitDefinition`` element."""

        # Check required attributes
        self._check_required_attributes(element, 'id')

        # Parse units and combine
        myokit_unit = myokit.units.dimensionless
        for unit in element.findall(self._path('listOfUnits', 'unit')):
            myokit_unit *= self._parse_unit(unit, model)

        # Store the unit
        sid = element.get('id')
        try:
            model.add_unit(sid, myokit_unit)
        except SBMLError as e:
            raise SBMLParsingError(
                'Unable to parse unit definition for "' + sid + '": ' + str(e),
                element)

    def _parse_unit(self, element, model):
        """Parses a ``unit`` element and returns a :class:`myokit.Unit`."""

        # Check required attributes
        self._check_required_attributes(element, 'kind')
        # The attributes 'multiplier', 'scale', and 'exponent' are required in
        # level 3, but are optional in level 2.

        # Get base unit (must be a predefined type)
        try:
            unit = model.base_unit(element.get('kind'))
        except SBMLError as e:
            raise SBMLParsingError(
                'Unable to parse unit kind: ' + str(e), element)

        # Parse remaining parts and return
        try:
            unit *= float(element.get('multiplier', default=1))
            unit *= 10 ** float(element.get('scale', default=0))
            unit **= float(element.get('exponent', default=1))
        except ValueError:
            raise SBMLParsingError('Unable to parse unit attributes', element)

        return unit

    def _path(self, *tags):
        """
        Returns a string created by prepending the namespace to each tag and
        adding forward slashes to separate.

        If a tag starts with a forward slash or period, no namespace will be
        prepended.
        """
        treated = []
        for tag in tags:
            if tag[:1] not in './':
                tag = '{' + self._ns + '}' + tag
            tag = tag.rstrip('/')
            treated.append(tag)
        return '/'.join(treated)

    def _tag(self, element):
        """
        Returns an element's name, but changes the syntax from ``{...}tag`` to
        ``sbml:tag`` for SBML elements.
        """
        ns, el = myokit.formats.xml.split(element.tag)

        # Replace known namespaces
        if ns is None:
            tag = el
        elif ns == self._ns:
            tag = 'sbml:' + el
        elif ns == NS_MATHML:
            tag = 'mathml:' + el
        else:
            tag = '{' + ns + '}' + el

        # Add id if known
        sid = element.get('id', element.get('unitsid'))
        if sid is not None:
            tag += '[@id=' + str(sid) + ']'

        return tag


class Quantity(object):
    """
    Base class for anything that has a numerical value in an SBML model, and
    can be set by rules, reactions, or initial assignments.
    """
    def __init__(self):

        # Expression for initial value (myokit.Expression)
        self._initial_value = None

        # Expression for value (myokit.Expression)
        self._value = None

        # Whether or not this is a rate
        self._is_rate = False

    def initial_value(self):
        """
        Returns an expression for this quantity's initial value, or ``None`` if
        not set.
        """
        return self._initial_value

    def rate(self):
        """
        Returns ``True`` only if this quantity's value is defined through a
        rate.
        """
        return self._is_rate

    def set_initial_value(self, value):
        """
        Sets a :class:`myokit.Expression` for this quantity's initial value.
        """
        self._initial_value = value

    def set_value(self, value, rate=False):
        """
        Sets a :class:`myokit.Expression` for this quantity's value.

        Arguments:

        ``value``
            An expression
        ``rate``
            Set to ``True`` if the expression gives the rate of change of this
            variable.

        """
        self._value = value
        self._is_rate = bool(rate)

    def value(self):
        """
        Returns an expression for this quantity's value, or ``None`` if not
        set.
        """
        return self._value


class Compartment(Quantity):
    """
    Represents a compartment in SBML; to create a compartment use
    :meth:`model.add_compartment()`.

    A compartment acts as a :class:`Quantity`, where the value represents the
    compartment's size.

    Arguments:

    ``model``
        The model this compartment is in.
    ``sid``
        This compartment's SId.

    """
    def __init__(self, model, sid):
        super(Compartment, self).__init__()
        self._model = model
        self._sid = sid

        self._spatial_dimensions = None
        self._size_units = None

    def set_spatial_dimensions(self, dimensions):
        """
        Sets the dimensionality of this compartment's size (usually 1, 2, or
        3); this is used to determine the units for the compartment size in the
        case that they are not explicitly set.
        """
        self._spatial_dimensions = float(dimensions)

    def set_size_units(self, units):
        """Sets the units for this compartment's size."""
        self._size_units = units

    def sid(self):
        """Returns this compartment's sid."""
        return self._sid

    def size_units(self):
        """
        Returns the units this compartment's size is in.

        If not set explicitly, the units will be retrieved from the model. If
        not set there either, units of dimensionless are returned.
        """
        if self._size_units is not None:
            return self._size_units

        if self._spatial_dimensions == 1:
            return self._model.length_units()
        elif self._spatial_dimensions == 2:
            return self._model.area_units()
        elif self._spatial_dimensions == 3:
            return self._model.volume_units()
        else:
            # Other values are allowed, but have no model attribute
            return myokit.units.dimensionless

    def spatial_dimensions(self):
        """
        Returns this compartment's spatial dimensions, or ``None`` if not set.
        """
        return self._spatial_dimensions

    def __str__(self):
        return '<Compartment ' + str(self._sid) + '>'


class Model(object):
    """
    Represents a model in SBML.

    Arguments:

    ``name``
        A user-friendly name for this model.

    """
    def __init__(self, name):

        self._name = name

        # Optional notes
        self._notes = None

        # Maps unit names to Unit objects
        # Units are the only things that can have a UnitSId
        self._units = {}

        # Used SIds
        self._sids = set()

        # Assignables: Compartments, species, species references, parameters
        self._assignables = {}

        # Compartments, parameters, species (all maps from sids to objects)
        self._compartments = collections.OrderedDict()
        self._parameters = collections.OrderedDict()
        self._species = collections.OrderedDict()

        # Reactions (map from sids to objects)
        self._reactions = collections.OrderedDict()

        # Default compartment size units
        self._length_units = myokit.units.dimensionless
        self._area_units = myokit.units.dimensionless
        self._volume_units = myokit.units.dimensionless

        # Default species amount units (not concentration)
        self._substance_units = myokit.units.dimensionless

        # Default time and extent units (kinetic law = extentUnits / timeUnits)
        self._time_units = myokit.units.dimensionless
        self._extent_units = myokit.units.dimensionless

        # Global conversion factor for species
        self._conversion_factor = None

    def add_compartment(self, sid):
        """Adds an SBML compartment to this model."""
        self._register_sid(sid)
        c = Compartment(self, sid)
        self._compartments[sid] = c
        self._assignables[sid] = c
        return c

    def add_parameter(self, sid):
        """Adds a parameter to this model."""
        self._register_sid(sid)
        p = Parameter(self, sid)
        self._parameters[sid] = p
        self._assignables[sid] = p
        return p

    def add_reaction(self, sid):
        """Adds a reaction to this model."""
        self._register_sid(sid)
        r = Reaction(self, sid)
        self._reactions[sid] = r
        return r

    def add_species(
            self, compartment, sid, amount=False, constant=False,
            boundary=False):
        """
        Adds a species to this model (located in the given ``compartment``).
        """
        self._register_sid(sid)
        s = Species(compartment, sid, amount, constant, boundary)
        self._species[sid] = s
        self._assignables[sid] = s

        return s

    def add_unit(self, unitsid, unit):
        """Adds a user unit with the given ``unitsid`` and myokit ``unit``."""
        if not _re_id.match(unitsid):
            raise SBMLError('Invalid UnitSId "' + str(unitsid) + '".')
        if unitsid in self._base_units or unitsid == 'celsius':
            raise SBMLError(
                'User unit overrides built-in unit: "' + str(unitsid) + '".')
        if unitsid in self._units:
            raise SBMLError(
                'Duplicate UnitSId: "' + str(unitsid) + '".')
        u = self._units[unitsid] = unit
        return u

    def area_units(self):
        """
        Returns the default compartment size units for 2-dimensional
        compartments, or dimensionless if not set.
        """
        return self._area_units

    def assignable(self, sid):
        """
        Returns the :class:`Compartment`, :class:`Species`,
        :class:`SpeciesReference`, or :class:`Parameter` referenced by the
        given SId.
        """
        return self._assignables[sid]

    def base_unit(self, unitsid):
        """
        Returns an SBML base unit, raises an :class:`SBMLError` if not
        supported.
        """
        # Check this base unit is supported
        if unitsid == 'celsius':
            raise SBMLError('The units "celsius" are not supported.')

        # Find and return
        return self._base_units[unitsid]

    def compartment(self, sid):
        """Returns the compartment with the given sid."""
        return self._compartments[sid]

    def conversion_factor(self):
        """
        Returns the :class:`Parameter` acting as global species conversion
        factor, or ``None`` if not set, see :meth:`Species.conversion_factor`.
        """
        return self._conversion_factor

    def extent_units(self):
        """
        Returns the default extent units in this model, or dimensionless if not
        set.
        """
        return self._extent_units

    def length_units(self):
        """
        Returns the default compartment size units for 1-dimensional
        compartments, or dimensionless if not set.
        """
        return self._length_units

    def myokit_model(self):
        """
        Converts this model to a :class:`myokit.Model`.

        SBML IDs are used for naming components and variables. If an ID starts
        with an underscore, the myokit name will be converter to
        `underscore<name>`.

        Compartments defined by the SBML file are mapped to myokit.Components.
        """

        # Create myokit model
        myokit_model = myokit.Model(self.name())

        # Add notes
        notes = self.notes()
        if notes:
            myokit_model.meta['desc'] = notes

        # Create reference container that links sid's to myokit objects
        component_references = {}
        variable_references = {}

        # Create additional container for species variables in amount
        # (Reactions define rates exclusively for amount of species)
        species_amount_references = {}

        # Create reference from Model expressions to myokit model expressions
        expression_references = {}

        # Instantiate component and variable objects first without assigning
        # RHS expressions. Myokit objects may have to be renamed, so
        # expressions are added in a second step.

        # Add SBML compartments to myokit model
        for sid, compartment in self._compartments.items():
            # Create component for compartment
            component = myokit_model.add_component_allow_renaming(
                convert_name(sid))

            # Add component to reference list
            component_references[sid] = component

            # Add compartment size to component
            var = component.add_variable_allow_renaming('size')

            # Set unit of size variable
            var.set_unit(compartment.size_units())

            # Add size variable to reference list
            variable_references[sid] = var
            expression_references[myokit.Name(compartment)] = myokit.Name(var)

        # Add myokit component to model to store time bound variable
        # and global parameters

        # Create default name of component
        myokit_component_name = 'myokit'

        # Add component
        component = myokit_model.add_component_allow_renaming(
            myokit_component_name)

        # Add myokit component to referece list
        component_references[myokit_component_name] = component

        # Add species to components
        for sid, species in self._species.items():
            # Get component from reference list
            compartment = species.compartment()
            compartment_sid = compartment.sid()
            component = component_references[compartment_sid]

            # Add species in amount to component
            # (needed for reactions, even if species is defined in
            # concentration)
            var = component.add_variable_allow_renaming(
                sid + '_amount')

            # Set unit of amount
            var.set_unit(species.substance_units())

            # Add reference to amount variable
            species_amount_references[sid] = var

            # Add species in concentration if measured in concentration
            if not species.amount():
                # Add species in concentration
                var = component.add_variable_allow_renaming(
                    sid + '_concentration')

                # Get myokit amount and size variable
                amount = species_amount_references[sid]
                size = variable_references[compartment_sid]

                # Set unit of concentration
                var.set_unit(amount.unit() / size.unit())

                # Define RHS of concentration as amount / size
                rhs = myokit.Divide(
                    myokit.Name(amount),
                    myokit.Name(size))
                var.set_rhs(rhs)

            # Add reference to species (either in amount or concentration)
            variable_references[sid] = var
            expression_references[myokit.Name(species)] = myokit.Name(var)

        # Add parameters to myokit component
        for sid, parameter in self._parameters.items():
            # Get myokit component
            component = component_references[myokit_component_name]

            # Add parameter to component
            var = component.add_variable_allow_renaming(sid)

            # Set unit of parameter
            var.set_unit(parameter.units())

            # Add reference to parameter
            variable_references[sid] = var
            expression_references[myokit.Name(parameter)] = myokit.Name(var)

        # Add stoichiometries from reactions
        for reaction in self._reactions.values():
            # Get all reactants and products (modifier do not have
            # stoichiometries)
            species_references = reaction.reactants() + reaction.products()

            for species_reference in species_references:
                sid = species_reference.sid()
                if sid is not None:
                    # Get component
                    species = species_reference.species()
                    compartment = species.compartment()
                    compartment_sid = compartment.sid()
                    component = component_references[compartment_sid]

                    # Add variable to component
                    var = component.add_variable_allow_renaming(sid)

                    # Set unit of variable
                    # (SBML defines stoichiometries as dimensionless)
                    var.set_unit(myokit.units.dimensionless)

                    # Add reference to variable
                    variable_references[sid] = var
                    expression_references[
                        myokit.Name(species_reference)] = myokit.Name(var)

        # Add time variable to myokit component
        component = component_references[myokit_component_name]
        var = component.add_variable_allow_renaming('time')

        # Bind time variable to time in myokit model
        var.set_binding('time')

        # Set time unit and initial value (SBML default: t=0)
        var.set_unit(self.time_units())
        var.set_rhs(0)

        # Add reference to time variable
        # (SBML referes to time by a csymbol:
        # 'http://www.sbml.org/sbml/symbols/time/')
        variable_references[
            'http://www.sbml.org/sbml/symbols/time'] = var

        # Set RHS of compartment sizes
        for sid, compartment in self._compartments.items():
            # Get myokit variable
            var = variable_references[sid]

            # Set initial value
            expr = compartment.initial_value()
            if expr is not None:
                var.set_rhs(expr.clone(
                    subst=expression_references))

            if compartment.rate():
                # Promote size to state variable
                var.promote(state_value=var.eval())

            # Set RHS
            # (assignmentRule overwrites initialAssignment)
            expr = compartment.value()
            if expr is not None:
                var.set_rhs(expr.clone(
                    subst=expression_references))

        # Set RHS of species (initalAssignemnt, assignmentRule, rateRule)
        # NOTE: Concentration is already set by (amount / size). So we
        # only set RHS of amount here.
        for sid, species in self._species.items():
            # Get myokit variable
            # We only adapt amount of species
            var = species_amount_references[sid]

            # Set initial value
            expr = species.initial_value()
            if expr is not None:
                # Need to convert initial value if
                # 1. the species is in amount but initial value units are not
                # 2. the species and the initial value is in concentration
                if species.amount() != species.units_initial_value():
                    # Get initial compartment size
                    compartment = species.compartment()
                    size = compartment.initial_value()

                    # Convert initial value from concentration to amount
                    expr = myokit.Multiply(expr, size)

                # Set initial value
                var.set_rhs(expr.clone(subst=expression_references))

            if species.rate():
                # Promote species to state variable
                var.promote(var.eval())

            # Set RHS (reactions are dealt with elsewhere)
            expr = species.value()
            if expr is not None:
                # Need to convert initial value if species is measured in
                # concentration (assignments match unit of measurement)
                if not species.amount():
                    # Get initial compartment size
                    compartment = species.compartment()
                    size = compartment.initial_value()

                    # Convert initial value from concentration to amount
                    expr = myokit.Multiply(expr, size)

                # Set initial value
                var.set_rhs(expr.clone(subst=expression_references))

        # Set RHS of parameters
        for sid, parameter in self._parameters.items():
            # Get myokit variable
            var = variable_references[sid]

            # Set initial value
            expr = parameter.initial_value()
            if expr is not None:
                var.set_rhs(expr.clone(
                    subst=expression_references))

            if parameter.rate():
                # Promote parameter to state variable
                var.promote(state_value=var.eval())

            # Set RHS
            # (assignmentRule overwrites initialAssignment)
            expr = parameter.value()
            if expr is not None:
                var.set_rhs(expr.clone(
                    subst=expression_references))

        # Set RHS of stoichiometries
        for reaction in self._reactions.values():
            # Get all reactants and products (modifier do not have
            # stoichiometries)
            species_references = reaction.reactants() + reaction.products()

            for species_reference in species_references:
                # Get sid
                sid = species_reference.sid()

                # Get stoichiometry variable
                try:
                    var = variable_references[sid]
                except KeyError:
                    continue

                # Set initial value
                expr = species_reference.initial_value()
                if expr is not None:
                    var.set_rhs(expr.clone(
                        subst=expression_references))

                if species_reference.rate():
                    # Promote stoichiometry to state variable
                    var.promote(state_value=var.eval())

                # Set RHS
                # (assignmentRule overwrites initialAssignment)
                expr = species_reference.value()
                if expr is not None:
                    var.set_rhs(expr.clone(
                        subst=expression_references))

        # Set RHS of species changed by reactions
        for reaction in self._reactions.values():
            # Get kinetic law
            kinetic_law = reaction.kinetic_law()

            if kinetic_law is None:
                # Skip to the next reaction
                continue

            for reactant in reaction.reactants():
                # Get species object
                species = reactant.species()

                if species.constant() or species.boundary():
                    # Species is not altered by the reaction
                    # Skip to the next reactant
                    continue

                # Instantiate rate expression
                expr = kinetic_law.clone()

                # Get stoichiometry of reactant
                try:
                    stoichiometry = myokit.Name(
                        variable_references[reactant.sid()])
                except KeyError:
                    stoichiometry = reactant.initial_value()

                if stoichiometry is not None:
                    # Weight rate expression by stoichiometry
                    expr = myokit.Multiply(stoichiometry, expr)

                factor = species.conversion_factor()
                if factor is not None:
                    # Get conversion factor variable
                    sid = factor.sid()
                    conversion_factor = variable_references[sid]

                    # Convert rate expression from units of reaction extent
                    # to amount units
                    conversion_factor = myokit.Name(conversion_factor)
                    expr = myokit.Multiply(conversion_factor, expr)

                # Get myokit amount variable
                var = species_amount_references[species.sid()]

                if not var.is_state():
                    # Promote amount to state variable
                    try:
                        var.promote(state_value=var.eval())
                        var.set_rhs(None)
                    except AttributeError:
                        var.promote()
                        var.set_rhs(None)

                if var.rhs().eval():
                    # Subtract rate contributions
                    # (Reaction removes species from compartment)
                    expr = myokit.Minus(var.rhs(), expr)
                else:
                    expr = myokit.PrefixMinus(expr)

                # Set RHS
                var.set_rhs(expr.clone(subst=expression_references))

            for product in reaction.products():
                # Get species object
                species = product.species()

                if species.constant() or species.boundary():
                    # Species is not altered by the reaction
                    # Skip to the next product
                    continue

                # Instantiate rate expression
                expr = kinetic_law.clone()

                # Get stoichiometry of product
                try:
                    stoichiometry = myokit.Name(
                        variable_references[product.sid()])
                except KeyError:
                    stoichiometry = product.initial_value()

                if stoichiometry is not None:
                    # Weight rate expression by stoichiometry
                    expr = myokit.Multiply(stoichiometry, expr)

                if species.conversion_factor():
                    # Get conversion factor variable
                    sid = species.conversion_factor().sid()
                    conversion_factor = variable_references[sid]

                    # Convert rate expression from units of reaction extent
                    # to amount units
                    conversion_factor = myokit.Name(conversion_factor)
                    expr = myokit.Multiply(conversion_factor, expr)

                # Get myokit amount variable
                var = species_amount_references[species.sid()]

                if not var.is_state():
                    # Promote amount to state variable
                    try:
                        var.promote(state_value=var.eval())
                        var.set_rhs(None)
                    except AttributeError:
                        var.promote()
                        var.set_rhs(None)

                if var.rhs().eval():
                    # Add rate contributions
                    expr = myokit.Plus(var.rhs(), expr)

                # Set RHS
                var.set_rhs(expr.clone(subst=expression_references))

        return myokit_model

    def name(self):
        """Returns this model's name."""
        return self._name

    def notes(self):
        """Returns a string of notes (if set), or ``None``."""
        return self._notes

    def parameter(self, sid):
        """Returns the :class:`Parameter` with the given id."""
        return self._parameters[sid]

    def reaction(self, sid):
        """Returns the :class:`Reaction` with the given id."""
        return self._reactions[sid]

    def _register_sid(self, sid):
        """
        Checks if the given SId is valid and available, registers it if it is,
        raises an error if it's not.
        """
        if not _re_id.match(sid):
            raise SBMLError('Invalid SId "' + str(sid) + '".')
        if sid in self._sids:
            raise SBMLError('Duplicate SId "' + str(sid) + '".')
        self._sids.add(sid)

    def set_area_units(self, units):
        """
        Sets the default compartment size units for 2-dimensional compartments.
        """
        self._area_units = units

    def set_conversion_factor(self, factor):
        """
        Sets a :class:`Parameter` as global conversion factor for species,
        see :meth:`Species.conversion_factor()`.
        """
        self._conversion_factor = factor

    def set_extent_units(self, units):
        """
        Sets the default units for "reaction extent", i.e. for the kinetic law
        equations in reactions.
        """
        self._extent_units = units

    def set_length_units(self, units):
        """
        Sets the default compartment size units for 1-dimensional compartments.
        """
        self._length_units = units

    def set_notes(self, notes=None):
        """Sets an optional notes string for this model."""
        self._notes = None if notes is None else str(notes)

    def set_substance_units(self, units):
        """Sets the default units for reaction amounts (not concentrations)."""
        self._substance_units = units

    def set_time_units(self, units):
        """Sets the time units used throughout the model."""
        self._time_units = units

    def set_volume_units(self, units):
        """
        Sets the default compartment size units for 3-dimensional compartments.
        """
        self._volume_units = units

    def __str__(self):
        if self._name is None:
            return '<SBMLModel>'
        return '<SBMLModel ' + str(self._name) + '>'

    def species(self, sid):
        """Returns the species with the given id."""
        return self._species[sid]

    def substance_units(self):
        """
        Returns the default units for reaction amounts (not concentrations), or
        dimensionless if not set.
        """
        return self._substance_units

    def time_units(self):
        """Returns the default units for time, or dimensionless if not set."""
        return self._time_units

    def unit(self, unitsid):
        """Returns a user-defined or predefined unit."""
        try:
            return self._units[unitsid]
        except KeyError:
            return self.base_unit(unitsid)

    def volume_units(self):
        """
        Returns the default compartment size units for 3-dimensional
        compartments, or dimensionless if not set.
        """
        return self._volume_units

    # SBML base units (except Celsius, because it's not defined in myokit)
    _base_units = {
        'ampere': myokit.units.A,
        'avogadro': myokit.parse_unit('1 (6.02214179e23)'),
        'becquerel': myokit.units.Bq,
        'candela': myokit.units.cd,
        'coulomb': myokit.units.C,
        'dimensionless': myokit.units.dimensionless,
        'farad': myokit.units.F,
        'gram': myokit.units.g,
        'gray': myokit.units.Gy,
        'henry': myokit.units.H,
        'hertz': myokit.units.Hz,
        'item': myokit.units.dimensionless,  # Myokit does not have item unit
        'joule': myokit.units.J,
        'katal': myokit.units.kat,
        'kelvin': myokit.units.K,
        'kilogram': myokit.units.kg,
        'liter': myokit.units.L,
        'litre': myokit.units.L,
        'lumen': myokit.units.lm,
        'lux': myokit.units.lux,
        'meter': myokit.units.m,
        'metre': myokit.units.m,
        'mole': myokit.units.mol,
        'newton': myokit.units.N,
        'ohm': myokit.units.ohm,
        'pascal': myokit.units.Pa,
        'radian': myokit.units.rad,
        'second': myokit.units.s,
        'siemens': myokit.units.S,
        'sievert': myokit.units.Sv,
        'steradian': myokit.units.sr,
        'tesla': myokit.units.T,
        'volt': myokit.units.V,
        'watt': myokit.units.W,
        'weber': myokit.units.Wb,
    }


class Parameter(Quantity):
    """
    Represents a parameter in SBML; to create a parameter use
    :meth:`model.add_parameter()`.

    Arguments:

    ``model``
        The model this parameter is in.
    ``sid``
        This parameter's SId.

    """
    def __init__(self, model, sid):
        super(Parameter, self).__init__()

        self._model = model
        self._sid = sid
        self._units = None

    def set_units(self, units):
        """Sets this parameters units to the given ``Units``."""
        self._units = units

    def sid(self):
        """Returns this parameter's sid."""
        return self._sid

    def __str__(self):
        return '<Parameter ' + str(self._sid) + '>'

    def units(self):
        """Returns the units this parameter is in, or ``None`` if not set."""
        return self._units


class Reaction(object):
    """
    Represents an SBML reaction; to create a reaction use
    :meth:`Model.add_reaction()`.

    Arguments:

    ``model``
        The model this reaction is in.
    ``sid``
        This reaction's SId.

    """
    def __init__(self, model, sid):
        self._model = model
        self._sid = sid

        # Reactants, reaction products, and modifiers: all as SpeciesReference
        self._reactants = []
        self._products = []
        self._modifiers = []

        # All species involved in this reaction (sid to object)
        self._species = collections.OrderedDict()

        # The kinetic law specifying this reaction's rate (if set)
        self._kinetic_law = None

    def add_modifier(self, species, sid=None):
        """Adds a modifier to this reaction and returns the created object."""
        if sid is not None:
            self._model._register_sid(sid)
        ref = ModifierSpeciesReference(species, sid)
        self._modifiers.append(ref)
        self._species[species.sid()] = species
        return ref

    def add_product(self, species, sid=None):
        """
        Adds a reaction product to this reaction and returns the created
        object.
        """
        if sid is not None:
            self._model._register_sid(sid)
        ref = SpeciesReference(species, sid)
        self._products.append(ref)
        self._species[species.sid()] = species
        if sid is not None:
            self._model._assignables[sid] = ref
        return ref

    def add_reactant(self, species, sid=None):
        """Adds a reactant to this reaction and returns the created object."""
        if sid is not None:
            self._model._register_sid(sid)
        ref = SpeciesReference(species, sid)
        self._reactants.append(ref)
        self._species[species.sid()] = species
        if sid is not None:
            self._model._assignables[sid] = ref
        return ref

    def kinetic_law(self):
        """
        Returns the kinetic law set for this reaction, or ``None`` if not set.
        """
        return self._kinetic_law

    def modifiers(self):
        """
        Returns the list of modifiers used by this reaction (as
        :class:`SpeciesReference` objects).
        """
        return list(self._modifiers)

    def products(self):
        """
        Returns the list of products created by this reaction (as
        :class:`SpeciesReference` objects).
        """
        return list(self._products)

    def reactants(self):
        """
        Returns the list of reactions used by this reaction (as
        :class:`SpeciesReference` objects).
        """
        return list(self._reactants)

    def set_kinetic_law(self, expression):
        """
        Sets this reaction's kinetic law (as a :class:`myokit.Expression`).
        """
        self._kinetic_law = expression

    def sid(self):
        """Returns this reaction's sid."""
        return self._sid

    def species(self, sid):
        """
        Finds and returns a :class:`Species` used in this reaction.
        """
        return self._species[sid]

    def __str__(self):
        return '<Reaction ' + str(self._sid) + '>'


class Species(Quantity):
    """
    Represents an SBML species; to create a species use
    :meth:`Compartment.add_species()`.

    Arguments:

    ``compartment``
        The :class:`Compartment` that this species is in.
    ``sid``
        This species's SId.
    ``amount``
        Whether this species value is represented as an amount (if false, it is
        represented as a concentration, which depends on the size of the
        compartment it is in).
    ``constant``
        Whether or not this species is constant.
    ``boundary``
        Whether or not this species is at the boundary of a reaction.

    """
    def __init__(self, compartment, sid, amount, constant, boundary):
        super(Species, self).__init__()

        self._compartment = compartment

        self._sid = sid
        self._amount = bool(amount)
        self._constant = bool(constant)
        self._boundary = bool(boundary)

        # Units for the amount, not the concentration
        self._units = None

        # Optional conversion factor from substance to extent units
        self._conversion_factor = None

        # Flag whether initial value is in correct units
        # (amount or concentration)
        self._units_initial_value = True

    def amount(self):
        """
        Returns ``True`` only if this species is defined as an amount (not a
        concentration).
        """
        return self._amount

    def boundary(self):
        """Returns ``True`` only if this species is at a reaction boundary."""
        return self._boundary

    def compartment(self):
        """Returns the :class:`Compartment` that this species is in."""
        return self._compartment

    def constant(self):
        """Returns ``True`` only if this species is constant."""
        return self._constant

    def conversion_factor(self):
        """
        Returns the :class:`Parameter` acting as conversion factor for this
        species, or the model default if not set (and ``None`` if that isn't
        set either).

        When calculating the rate of change of the species, multiplying by this
        factor converts from the species units to "number of reaction events".
        """
        if self._conversion_factor is not None:
            return self._conversion_factor
        return self._compartment._model.conversion_factor()

    def set_conversion_factor(self, factor):
        """
        Sets a :class:`Parameter` as conversion factor for this species,
        see :meth:`conversion_factor()`.
        """
        self._conversion_factor = factor

    def set_substance_units(self, units):
        """Sets the units this species amount (not concentration) is in."""
        self._units = units

    def set_units_initial_value(self, flag):
        """Sets a flag whether the units of the initial value are correct."""
        self._units_initial_value = flag

    def units_initial_value(self):
        """
        Returns a boolean flag whether the initial value has the correct units.
        """
        return self._units_initial_value

    def sid(self):
        """Returns this species's sid."""
        return self._sid

    def __str__(self):
        return '<Species ' + str(self._sid) + '>'

    def substance_units(self):
        """
        Returns the units an amount of this species (not a concentration) is
        in, or the model default if not set.
        """
        if self._units is not None:
            return self._units
        return self._compartment._model.substance_units()


class SpeciesReference(Quantity):
    """
    Represents a reference to a reactant or product species in an SBML
    reaction.

    A species reference acts as a :class:`Quantity`, where the value represents
    the species' stoichiometry.
    """
    REACTANT = 0
    PRODUCT = 1
    MODIFIER = 2

    def __init__(self, species, sid=None):
        super(SpeciesReference, self).__init__()

        self._species = species
        self._sid = sid

    def species(self):
        """Returns the species this object refers to."""
        return self._species

    def sid(self):
        """Returns this species reference's SId, or ``None`` if not set."""
        return self._sid


class ModifierSpeciesReference(object):
    """Represents a reference to a modifier species in an SBML reaction."""

    def __init__(self, species, sid=None):
        super(ModifierSpeciesReference, self).__init__()

        self._species = species
        self._sid = sid

    def species(self):
        """Returns the species this object refers to."""
        return self._species

    def sid(self):
        """Returns this species reference's SId, or ``None`` if not set."""
        return self._sid


def convert_name(name):
    """
    Converts SBML name to myokit compatible name.

    Unlike the SBML SId type, Myokit names cannot start with an
    underscore.
    """
    if name[:1] == '_':
        name = 'underscore_' + name
    return name

