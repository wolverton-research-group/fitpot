from utils import *
import json
import os
import random
import subprocess
import tempfile
import random
from numpy import linspace

class Genome:
    cont_vars = {}
    functions = {}
    disc_vars = {}
    constraints = []
    constants = {}

    def __init__(self):
        self.genes = {}
        self.randomize()
        self.fitness = 1e9
        while not self.valid:
            self.randomize()

    @property
    def generator_string(self):
        raise NotImplementedError

    @property
    def values(self):
        values = {}
        for k in self.cont_vars.keys():
            values[k] = self.genes[k] 
        for k, func in self.functions.items():
            values[k] = func(self)
        for k, v in self.disc_vars.items():
            if self.genes[k]:
                values[k] = v[0]
            else:
                values[k] = v[1]
        return values

    @property
    def valid(self):
        keys = self.cont_vars.keys() 
        keys += self.functions.keys()
        keys += self.disc_vars.keys()
        keys += self.constants.keys()
        if not len(set(keys)) == len(keys):
            return False
        return all( constraint(self) for constraint in self.constraints )

    def randomize(self):
        for k, vals in self.cont_vars.items():
            self.genes[k] = random.choice(
                    linspace(self.cont_vars[k][0], self.cont_vars[k][1], 1e3))
        for k in self.disc_vars:
            self.genes[k] = random.choice(self.disc_vars[k])

    def serial_repr(self):
        return json.dumps({
            'functions':self.functions,
            'cont_vars':self.cont_vars,
            'constants':self.constants,
            'cont_vars':self.cont_vars,
            'constraints':self.constraints,
            'genes':self.cont_vars})

    def load(self, jstring):
        jdict = json.loads(jstring)
        self.functions = jdict['functions']
        self.cont_vars = jdict['cont_vars']
        self.constraints = jdict['constraints']
        self.constants = jdict['constants']
        self.genes = jdict['genes']

    ### IPR specific parts
    @property
    def pot_string(self):
        return self.template.format(
                **self.values)

    def gulp(self, data):
        gulp_instr = 'conp gradient\n'
        gulp_instr += 'title\nIPR generated gulp input\nend\n'
        gulp_instr += 'vectors\n'
        for line in data.cell:
            gulp_instr += '%0.8f %0.8f %0.8f\n' % tuple(line)
        gulp_instr += '0 0 0 0 0 0\n'
        gulp_instr += 'fractional\n'
        for atom in data.coords:
            gulp_instr += '%s core %0.8f %0.8f %0.8f 1 0 0 0\n' % (
                    (atom[0],) + tuple(atom[1]))
            if 'shellmode' in self.constants:
                gulp_instr += '%s shell %0.8f %0.8f %0.8f 1 1 1 1\n' % (
                    (atom[0],) + tuple(atom[1]))
        gulp_instr += self.pot_string

        gulp_list = gulp_instr.split('\n')
        safe_gulp_list = []
        for line in gulp_list:
            safe_gulp_list.append(ensure_length(line))
        return '\n'.join(safe_gulp_list)

    def save(self, filename):
        f = open(filename,'w')
        f.write(self.pot_string)
        f.close()
