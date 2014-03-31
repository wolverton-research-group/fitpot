from genome import Genome
from numpy import linspace

__doc__='''
Fields are part of the basis of genome factories.

Field: 
    -name           :   some desription, i.e. many_body
    -container      :   'manybody'
    -cont_vars      :   {'key': [min_val, max_val] ... }
    -disc_vars      :   {'key': list_of_values ... }
    -functions      :   {'key': function_generator .. }
    -constraint     :   [ constraint_generators ... ]

    Generators within a field level define relationships between genes within a
    single line. i.e. r_cutoff for lennard jones is 2.5*sigma. Generators must
    take the field's prefix as an argument, and return a function that accepts
    an organism as its only input, and returns a value as its only output.
'''

class Field(object):
    _container = ''
    _field_prefix = ''
    _species_len = 0
    _template = ''

    _constants = {}
    _cont_vars = {}
    _functions = {}
    _disc_vars = {}
    _constraints = []

    def __init__(self, species, **kwargs):
        self.species = species
        self.options = kwargs
        assert self._species_len == len(species)

    @property
    def prefix(self):
        elts = [ spec.split()[0] for spec in self.species ]
        return self._field_prefix+'_'+'_'.join( e.lower() for e in elts )

    @property
    def container(self):
        return self._container

    @property
    def constants(self):
        return self._constants

    @property
    def cont_vars(self):
        vars = {}
        for k, v in self._cont_vars.items():
            vars[self.prefix+'_'+k] = v
        return vars

    @property
    def functions(self):
        funcs = {}
        for k, v in self._functions.items():
            funcs[self.prefix+'_'+k] = v(self.prefix)
        return funcs

    @property
    def disc_vars(self):
        togs = {}
        for k, v in self._disc_vars.items():
            togs[self.prefix+'_'+k] = v
        return togs

    @property
    def constraints(self):
        constrs = []
        for c in self._constraints:
            constrs.append(c(self.prefix))
        return constrs

    @property
    def keys(self):
        keys = []
        keys += self._cont_vars.keys()
        keys += self._functions.keys()
        keys += self._disc_vars.keys()
        keys += self._constants.keys()
        assert len(keys) == len(set(keys))
        return dict( (k,'{'+self.prefix+'_'+k+'}') for k in keys )

    @property
    def template(self):
        return ' '.join(self.species)+' '+self._template.format(**self.keys)

    def condition(self):
        return

class many_body(Field):
    container = 'manybody'
    _template='0.0 12.0'

class eam_functional_power_6(Field):
    container = 'eam_functional power 6'
    _field_prefix = 'efp6'
    _species_len = 1

    _template = '{A_1}'
    _cont_vars = {'A_1':[1e2, 5e3]}

class eam_density_sqrt(Field):
    container = 'eam_density square_root'
    _field_prefix = 'edsr'
    _species_len = 1

    _template = '{C}'
    _cont_vars=  {'C':[0, 1e6]}

class eam_repulsive(Field):
    container = 'lennard 12 6'
    _field_prefix = 'eam_2br'
    _species_len = 2

    _template = '{repel} 0.0 0.0 12.0'
    _cont_vars = {'repel':[1e2, 1e4]}


class lennard_jones(Field):
    container = 'lennard 12 6 epsilon'

    def scale_factory(prefix):
        def scale(pot):
            return 2.5*pot.genes[prefix+'_sig']
        return scale

    _field_prefix = 'lj'
    _species_len = 2

    _template = '{eps} {sig} 0.0 {scaled}'
    #_template = '{eps} {sig} 0.0 6.0'
    _functions = {'scaled':scale_factory}
    _cont_vars = {'eps':[0,100],
            'sig':[1.0,2.0]}
