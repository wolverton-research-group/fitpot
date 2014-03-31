from ipr import *
from ipr.library import LennardJones#, SuttonChen
from ipr.fields import many_body, eam_functional_power_6, eam_density_sqrt

### many flexible and powerful ways to create genomes
#genome = LennardJones(elements=['Al']).genome

#genome = LennardJones(elements=['Al','Ni'], 
#        cross_terms='geometric_mean').genome

run = Optimizer()
run.load_data('/home/sjk648/tmp/Al')

lj = LennardJones(elements=['Al'])

run.genome = lj.genome

#run()

#genome = genome_factory(
#        many_body(elt1='Al', elt2='Al'),
#        many_body(elt1='Al', elt2='Ni'),
#        many_body(elt1='Ni', elt2='Ni'),
#        eam_functional_power_6(elt='Al'),
#        eam_density_sqrt(elt='Al'),
#        eam_functional_power_6(elt='Ni'),
#        eam_density_sqrt(elt='Ni'))

#genome = SuttonChen(elt1='Al', elt2='Ni')
