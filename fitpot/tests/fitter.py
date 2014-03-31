from ipr import *
from ipr.library import LennardJones

run = Optimizer()
run.load_data('absolute path to data')

LJ = LennardJones(elements=['H','Al','Na'])

def charge_balance(pot):
    return -(7*pot.genes['H_core']+2*pot.genes['Al_core'])

LJ.species['H core'] = [-1, 0]
LJ.species['Al core'] = [0, 3]
LJ.species['Na core'] = charge_balance

run.genome = LJ.genome

run()
