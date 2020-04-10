#
# Imports a model from an SBML file.
# Only partial SBML support (Based on SBML level 3) is provided.
#
# This file is part of Myokit.
# See http://myokit.org for copyright, sharing, and licensing details.
#
from __future__ import absolute_import, division
from __future__ import print_function, unicode_literals

import xml.etree.ElementTree as ET
import xml.dom.minidom
import os
import re

import myokit
import myokit.units
import myokit.formats
from myokit.mxml import html2ascii
from myokit.mxml import dom_child, dom_next
from myokit.formats.mathml import parse_mathml_dom, parse_mathml_etree


info = """
Loads a Model definition from an SBML file. Warning: This importer hasn't been
extensively tested.
"""


class SBMLImporter(myokit.formats.Importer):
    """
    This:class:`Importer <myokit.formats.Importer>` load model definitions
    from files in SBML format.
    """
    def __init__(self):
        super(SBMLImporter, self).__init__()
        self.re_name = re.compile(r'^[a-zA-Z]+[a-zA-Z0-9_]*$')
        self.re_alpha = re.compile(r'[\W]+')
        self.re_white = re.compile(r'[ \t\f\n\r]+')
        self.units = {}

    def _convert_name(self, name):
        """
        Converts a name to something acceptable to myokit.
        """
        if not self.re_name.match(name):
            org_name = name
            name = self.re_white.sub('_', name)
            name = self.re_alpha.sub('_', name)
            if not self.re_name.match(name):
                name = 'x_' + name
            self.logger().warn(
                'Converting name <' + org_name + '> to <' + name + '>.')
        return name

    def _convert_unit(self, unit):
        """
        Converts an SBML unit to a myokit one using the lookup tables generated
        when parsing the XML file.
        """
        if unit in self.units:
            return self.units[unit]
        elif unit in unit_map:
            return unit_map[unit]
        else:      # pragma: no cover
            raise SBMLError('Unit not recognized: ' + str(unit))

    def info(self):
        return info

    def model(self, path):
        """
        Returns myokit model specified by the SBML fil provided.
        """
        # Get logger
        log = self.logger()

        # Parse xml file
        path = os.path.abspath(os.path.expanduser(path))
        tree = ET.parse(path)
        root = tree.getroot()

        # get SBML namespace
        sbml_version = self._get_sbml_version(root)
        self.ns = self._get_namespaces(sbml_version, log)

        # Get model node
        xmodel = root[0]
        if xmodel.get('name'):
            name = str(xmodel.get('name'))
        elif xmodel.get('id'):
            name = str(xmodel.get('id'))
        else:
            name = 'Imported SBML model'

        # Create myokit model
        model = myokit.Model(self._convert_name(name))
        log.log('Reading model "' + model.meta['name'] + '"')

        # Create one giant component to hold all variables
        comp = model.add_component('sbml')

        # Create table of variable names
        refs = {}

        # Handle notes, if given
        x = xmodel.find(self.ns['sbml'] + 'notes')
        if x:
            log.log('Converting <model> notes to ascii')
            model.meta['desc'] = html2ascii(ET.tostring(x, encoding='unicode'),
                                            width=75
                                            )
            # width = 79 - 4 for tab!

        # Warn about missing functionality
        x = xmodel.find(self.ns['sbml'] + 'listOfCompartments')
        if x:   # pragma: no cover
            log.warn('Compartments are not supported.')
        x = xmodel.find(self.ns['sbml'] + 'listOfSpecies')
        if x:   # pragma: no cover
            log.warn('Species are not supported.')
        x = xmodel.find(self.ns['sbml'] + 'listOfConstraints')
        if x:   # pragma: no cover
            log.warn('Constraints are not supported.')
        x = xmodel.find(self.ns['sbml'] + 'listOfReactions')
        if x:   # pragma: no cover
            log.warn('Reactions are not supported.')
        x = xmodel.find(self.ns['sbml'] + 'listOfEvents')
        if x:   # pragma: no cover
            log.warn('Events are not supported.')

        # Ignore custom functions
        x = xmodel.find(self.ns['sbml'] + 'listOfFunctionDefinitions')
        if x:   # pragma: no cover
            log.warn('Custom math functions are not (yet) implemented.')

        # Parse custom units
        x = xmodel.find(self.ns['sbml'] + 'listOfUnitDefinitions')
        if x:
            self._parse_units(model, comp, x, self.ns['sbml'])

        # Parse parameters (constants + parameters)
        x = xmodel.find(self.ns['sbml'] + 'listOfParameters')
        if x:
            self._parse_parameters(model, comp, refs, x, self.ns['sbml'])

        # Parse rules (equations)
        x = xmodel.find(self.ns['sbml'] + 'listOfRules')
        if x:
            self._parse_rules(model, comp, refs, x, self.ns['sbml'])

        # Parse extra initial assignments
        x = xmodel.find(self.ns['sbml'] + 'listOfInitialAssignments')
        if x:
            self._parse_initial_assignments(model,
                                            comp,
                                            refs,
                                            x,
                                            self.ns['sbml']
                                            )

        # Write warnings to log
        log.log_warnings()

        # Run model validation, order variables etc
        try:
            model.validate()
        except myokit.IntegrityError as e:
            log.log_line()
            log.log('WARNING: Integrity error found in model:')
            log.log(str(e))
            log.log_line()

        # Return finished model
        return model

    def _get_namespaces(self, sbml_version: str, log):
        """
        Defines namespaces for supported SBML versions. Supported versions are:
        level 2, version 3, 4, 5.
        """
        ns = dict()
        if sbml_version == "{http://www.sbml.org/sbml/level2/version3}":
            # SBML
            ns['sbml'] = "{http://www.sbml.org/sbml/level2/version3}"
            # MathML
            ns['mathml'] = "{http://www.w3.org/1998/Math/MathML}"
        elif sbml_version == "{http://www.sbml.org/sbml/level2/version4}":
            # SBML
            ns['sbml'] = "{http://www.sbml.org/sbml/level2/version4}"
            # MathML
            ns['mathml'] = "{http://www.w3.org/1998/Math/MathML}"
        elif sbml_version == "{http://www.sbml.org/sbml/level2/version5}":
            # SBML
            ns['sbml'] = "{http://www.sbml.org/sbml/level2/version5}"
            # MathML
            ns['mathml'] = "{http://www.w3.org/1998/Math/MathML}"
        elif sbml_version == "{http://www.sbml.org/sbml/level3/version1}":
            # SBML
            ns['sbml'] = "{http://www.sbml.org/sbml/level3/version1}"
            # MathML
            ns['mathml'] = "{http://www.w3.org/1998/Math/MathML}"
        elif sbml_version == "{http://www.sbml.org/sbml/level3/version2}":
            # SBML
            ns['sbml'] = "{http://www.sbml.org/sbml/level3/version2}"
            # MathML
            ns['mathml'] = "{http://www.w3.org/1998/Math/MathML}"
        else:
            # SBML
            ns['sbml'] = sbml_version
            # MathML
            ns['mathml'] = "{http://www.w3.org/1998/Math/MathML}"

            # log import warning
            log.warn('The SBML version %s has not been tested. The model'
                     % sbml_version + ' may not be imported correctly.'
                     )

        return ns

    def _get_sbml_version(self, root):
        """
        Returns the SBML version of the file.
        """
        m = re.match(r'\{.*\}', root.tag)

        return m.group(0) if m else ''

    def _parse_initial_assignments(self, model, comp, refs, node, ns):
        """
        Parses any initial values specified outside of the rules section.
        """
        nodes = node.findall(ns + 'initialAssignment')
        # get mathml ns
        mathml_ns = self.ns['mathml']
        # iterate through initial assignments
        for node in nodes:
            var = str(node.get('symbol')).strip()
            var = self._convert_name(var)
            if var in comp:
                self.logger().log(
                    'Parsing initial assignment for "' + var + '".')
                var = comp[var]
                # get child
                child = node.find(mathml_ns + 'math')
                if child:
                    expr = parse_mathml_etree(
                        child,
                        lambda x, y: myokit.Name(refs[x]),
                        lambda x, y: myokit.Number(x)
                    )

                    if var.is_state():
                        # Initial value
                        var.set_state_value(expr)
                    else:
                        # Change of value
                        var.set_rhs(expr)
            else:   # pragma: no cover
                raise SBMLError(   # pragma: no cover
                    'Initial assignment found for unknown parameter <' + var
                    + '>.')

    def _parse_parameters(self, model, comp, refs, node, ns):
        """
        Parses parameters
        """
        nodes = node.findall(ns + 'parameter')
        for node in nodes:
            # Create variable
            org_name = str(node.get('id'))
            name = self._convert_name(org_name)
            self.logger().log('Found parameter "' + name + '"')
            if name in comp:    # pragma: no cover
                self.logger().warn(
                    'Skipping duplicate parameter name: ' + str(name))
            else:
                # Create variable
                unit = None
                if node.get('units'):
                    foreign_unit = node.get('units')
                    if foreign_unit:
                        unit = self._convert_unit(foreign_unit)
                value = None
                if node.get('value'):
                    value = node.get('value')
                var = comp.add_variable(name)
                var.set_unit(unit)
                var.set_rhs(value)

                # Store reference to variable
                refs[org_name] = refs[name] = var

    def _parse_rules(self, model, comp, refs, node, ns):
        """
        Parses the rules (equations) in this model
        """
        parent = node
        # Create variables with assignment rules (all except derivatives)
        nodes = parent.findall(ns + 'assignmentRule')
        # get MathML ns
        mathml_ns = self.ns['mathml']
        # iterate through assignment rules
        for node in nodes:
            var = self._convert_name(
                str(node.get('variable')).strip())
            if var in comp:
                self.logger().log(
                    'Parsing assignment rule for <' + str(var) + '>.')
                var = comp[var]
                # get child
                child = node.find(mathml_ns + 'math')
                # add expression to model
                if child:
                    var.set_rhs(parse_mathml_etree(
                        child,
                        lambda x, y: myokit.Name(refs[x]),
                        lambda x, y: myokit.Number(x)
                    )
                    )
            else:
                raise SBMLError(   # pragma: no cover
                    'Assignment found for unknown parameter: "' + var + '".')

        # Create variables with rate rules (states)
        nodes = parent.findall(ns + 'rateRule')
        for node in nodes:
            var = self._convert_name(
                str(node.get('variable')).strip())
            if var in comp:
                self.logger().log('Parsing rate rule for <' + var + '>.')
                var = comp[var]
                ini = var.rhs()
                ini = ini.eval() if ini else 0
                var.promote(ini)
                # get child
                child = node.find(mathml_ns + 'math')
                # add expression to model
                if child:
                    var.set_rhs(parse_mathml_etree(
                        child,
                        lambda x, y: myokit.Name(refs[x]),
                        lambda x, y: myokit.Number(x)
                    )
                    )
            else:
                raise SBMLError(   # pragma: no cover
                    'Derivative found for unknown parameter: <' + var + '>.')

    def _parse_units(self, model, comp, node, ns):
        """
        Parses custom unit definitions, creating a look-up table that can be
        used to convert these units to myokit ones.
        """
        nodes = node.findall(ns + 'unitDefinition')
        for node in nodes:
            name = node.get('id')
            self.logger().log('Parsing unit definition for "' + name + '".')
            unit = myokit.units.dimensionless
            units = node.find(ns + 'listOfUnits')
            units = units.findall(ns + 'unit')
            for node2 in units:
                kind = str(node2.get('kind')).strip()
                u2 = self._convert_unit(kind)
                if node2.get('multiplier'):
                    m = float(node2.get('multiplier'))
                else:
                    m = 1.0
                if node2.get('scale'):
                    m *= 10 ** float(node2.get('scale'))
                u2 *= m
                if node2.get('exponent'):
                    u2 **= float(node2.get('exponent'))
                unit *= u2
            self.units[name] = unit

    def supports_model(self):
        return True


class SBMLError(myokit.ImportError):
    """
    Thrown if an error occurs when importing SBML
    """


unit_map = {
    'dimensionless': myokit.units.dimensionless,
    'ampere': myokit.units.A,
    'farad': myokit.units.F,
    'katal': myokit.units.kat,
    'lux': myokit.units.lux,
    'pascal': myokit.units.Pa,
    'tesla': myokit.units.T,
    'becquerel': myokit.units.Bq,
    'gram': myokit.units.g,
    'kelvin': myokit.units.K,
    'meter': myokit.units.m,
    'radian': myokit.units.rad,
    'volt': myokit.units.V,
    'candela': myokit.units.cd,
    'gray': myokit.units.Gy,
    'kilogram': myokit.units.kg,
    'metre': myokit.units.m,
    'second': myokit.units.s,
    'watt': myokit.units.W,
    'celsius': myokit.units.C,
    'henry': myokit.units.H,
    'liter': myokit.units.L,
    'mole': myokit.units.mol,
    'siemens': myokit.units.S,
    'weber': myokit.units.Wb,
    'coulomb': myokit.units.C,
    'hertz': myokit.units.Hz,
    'litre': myokit.units.L,
    'newton': myokit.units.N,
    'sievert': myokit.units.Sv,
    'joule': myokit.units.J,
    'lumen': myokit.units.lm,
    'ohm': myokit.units.R,
    'steradian': myokit.units.sr,
}
