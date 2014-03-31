from ipr import *
from numpy import linspace

run = Optimizer()

run.load_data('/home/sjk648/tmp/Al')

genome = Genome

genome.template = '''
species
Al core 0.0

{eam}manybody
{eam}Al core Al core 0 12
{lennard}lennard 6 12
{lennard}Al core Al core {eps} {sigma} 0.0 {scaled}

{eam}eam_density power {power}
{eam}Al core {C}
{eam}
{eam}eam_functional square_root
{eam}Al core {A_1}
'''

### Variables are given as a name, and a list of values it can take
genome.cont_vars = {
        'eps': [-5,0], 
        'sigma':[1,6],
        'C':[0,1e6],
        'A_1':[-5,0],
        'power_var':[1,10]}

### Toggles are given as a name, and a list of length 2, with the first element
###     being the value if True, the second the value if False
genome.disc_vars = {
        'lennard':['#',''],
        'eam':['#','']}


### Functions can be predifined and plugged in, or given as lambda functions. 

## explicitly
def scaled(pot):
    return 2.5*pot.genes['sigma']

## or with lambda functions
genome.functions = {
        'scaled':scaled,
        'power':lambda x: int(round(x.genes['power_var']))}


### Constraints, are basically functions, but are evaluated at organism
###     creation to ensure, for example, r_0 < r_1, or two types of potentials
###     aren't turned on at onces
genome.constraints = [
        lambda x: x.genes['lennard'] != x.genes['eam']]

run.genome = genome

run()
