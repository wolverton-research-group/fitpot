from genome import Genome
from fields import *
import itertools
from collections import defaultdict

__doc__='''
Base class for genome factories.
FunctionalForm:
    -name = crude description

    -fields         :   [ fields ... ]
    -cont_vars      :   {'key': list_of_values ... }
    -functions      :   {'key': function_generator ... }
    -disc_vars        :   {'key' : disc_vars ...} 
    -constraints    :   [ constraint_generators ... ]

    Generators at the factory level supercede field level definitions, and  
'''


class GenomeFactory(object):
    name = None
    fields = []
    cont_vars = {}
    functions = {}
    disc_vars = {}
    constants = {}
    constraints = []

    @property
    def genome(self):
        genome = Genome
        for field in self.fields:
            self.containers[field.container].append(field)            

        genome.template = 'species\n'
        for spec, charge in self.species.items():
            if isinstance(charge, list):
                genome.cont_vars[spec.replace(' ','_')] = charge
                genome.template += '%s {%s}\n' % (spec, spec.replace(' ','_'))
            elif type(charge) == type(genome_factory):
                genome.functions[spec.replace(' ','_')] = charge
                genome.template += '%s {%s}\n' % (spec, spec.replace(' ','_'))
            else:
                genome.template += '%s %s\n' % (spec, charge)

        for cont, fields in self.containers.items():
            genome.template += '%s\n' % cont
            for f in fields:
                genome.template += '     %s\n' % f.template
                genome.cont_vars.update(f.cont_vars)
                genome.functions.update(f.functions)
                genome.disc_vars.update(f.disc_vars)
                genome.constants.update(f.constants)
                genome.constraints += f.constraints
        return genome

def genome_factory(cont_vars={},
        functions={},
        disc_vars = {},
        constants = {},
        constraints=[], *args):
    gf = GenomeFactory
    gf.fields = args
    gf.functions = functions
    gf.disc_vars = disc_vars
    gf.constraints = constraints
    gf.constants = constants
    return gf.genome

class LennardJones(GenomeFactory):
    name = 'Lennard Jones'
    fields = []
    species = {}
    containers = defaultdict(list)

    _cont_vars = {}
    _functions = {}
    _disc_vars = {}
    _constraints = []

    def __init__(self, elements=[], cross_terms='fit'):
        for elt in elements:
            self.fields.append(lennard_jones([elt+' core', elt+' core']))
            self.species['%s core' % elt] = 0
        if cross_terms == 'fit':
            for e1, e2 in itertools.combinations(elements, 2):
                self.fields.append(lennard_jones([e1+' core', e2+' core']))
