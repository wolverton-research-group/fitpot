from ipr import *
from numpy import linspace

run = Optimizer()

run.load_data('/home/sjk648/tmp/Al')

genome = Genome

genome.template = '''
species
Al core 0.0

manybody
Al core Al core 0.0 12.0

eam_density power 6
Al core {A_1}

eam_functional square_root
Al core 1.0

## repulive two body components
lennard 7 6
Al core Al core {repulse} 0.0 0.0 12.0
'''

### Variables are given as a name, and a list of values it can take
genome.cont_vars = {
        'repulse':[1e2, 1e3],
        'A_1':[5e2, 5e3]}

### Toggles are given as a name, and a list of length 2, with the first element
###     being the value if True, the second the value if False
genome.disc_vars = {}

### Functions can be predifined and plugged in, or given as lambda functions. 

## or with lambda functions
genome.functions = {}

### Constraints, are basically functions, but are evaluated at organism
###     creation to ensure, for example, r_0 < r_1, or two types of potentials
###     aren't turned on at onces
genome.constraints = []
run.genome = genome

run.first_generation_factor = 10
run.fit_size = 5
run.p_mutate = 0.25
run.pop_size = 100
run.tourn_size = 5

run()
