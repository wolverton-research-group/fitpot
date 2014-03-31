# ipr.model

import random
import os
import sys
import time
import tempfile
import itertools
import multiprocessing as mp
import subprocess
import numpy as np
import csv
from utils import *
from config import *
from data import *
from collections import defaultdict
#from analysis import *

def gulp_call(bundle):
    org_ind, data_ind, gulp_instr = bundle
    tmp = tempfile.mkdtemp(dir='/dev/shm')
    cwd = os.getcwd()
    os.chdir(tmp)
    target = tempfile.mktemp(dir=tmp, suffix='.frc')
    gulp_instr += 'output frc '+target
    p1 = subprocess.Popen(['echo',gulp_instr],
            stdout=subprocess.PIPE)
    p2 = subprocess.Popen(GULP_CMD,
            #stdout=subprocess.PIPE,
            stdout=open('/dev/null','w'),
            stderr=subprocess.PIPE,
            stdin=p1.stdout)
    p2.wait()
    os.chdir(cwd)
    out, err = p2.communicate()
    try:
        frcout = open(target,'r')
        result = read_frcout(frcout)
        frcout.close()
        os.remove(target)
        os.rmdir(tmp)
    except:
        return (org_ind, data_ind, False)
    return (org_ind, data_ind, result)


class Optimizer:
    '''
    Optimizer objects hold all of the information and instructions for fitting
    the potential.

    [Fit Options]
    -fit_size[default=20]: how many data points to fit to
    -test_size[default=0]: how many data points to test against independently
    -pop_size[default=100]: how many organisms per generation
    -weights[default={'energy':10, 'stress':1, 'force':1}]: a dictionary
    containing the relative weights of force, stress and energy errors

    [GA Options]
    -tourn_size[default=5]: how many organisms to compete for procreation
    -f_replace[default=1.0]: how many organisms to keep back from the previous
    generation, to 'seed' the next one
    -p_mutate[default=0.1]: chance of mutating any given gene during breeding.
    -first_generation_size[default=500]: size of the first generation to get a
    pool of reasonably good structures to start with.

    [Extras]
    -rescale_threshold[default=
    -sliding_threshold[default=0.25]
    '''

    fit_size = 20
    tourn_size = 5
    pop_size = 100
    f_replace = 1.0
    p_mutate = 0.1
    first_generation_size = 500
    weights = {'energy':10.0,
            'stress':1.0,
            'force':1.0}
    rescale_threshold = 0.05
    sliding_threshold = 0.25


    def __init__(self, genome=None):
        ### core attributes
        self.genome = genome

        ### run cont_vars
        self.data = {}  ## all available data
        self.fit_set = [] ## data to fit to
        self.test_set = [] ## independent test set
        self.organisms = {}
        self.generations = []
        self.results = []

    def load_data(self, path):
        data = Data.read_path(path)
        for d in data:
            d.id = len(self.data)
            self.data[d.id] = d

    ### Evaluate fitness

    @property
    def best(self):
        evaled = [ self.organisms[o] for o in self.generations[-1]]
        return sorted( evaled, key=lambda x: x.fitness )[0]

    @property
    def genes(self):
        genes = []
        for org in self.organisms.values():
            genes.append(org[k] for k in self.keys)
        genes = np.array(genes)
        return genes.T

    @property
    def keys(self):
        if not self.organisms:
            self.create_organism()
        return self.organisms[0].genes.keys()

    def evaluate(self, org, data):
        if isinstance(org, int):
            org_ind = org
            org = self.organisms[org]
        else:
            org_ind = org.id
        if isinstance(data, int):
            data_ind = data
            data = self.data[data_ind]
        else:
            data_ind = data.id

        instr = org.gulp(data)
        return gulp_call((org_ind, data_ind, instr))

    def bulk_evaluate(self, 
            datas=None, 
            orgs=None,
            verbosity=0):

        if not datas:
            datas = list(self.fit_set)
        if not orgs:
            orgs = list(self.generations[-1])

        p = mp.Pool()
        tests =  itertools.product(orgs, datas)
        todo = [ (org, data, self.organisms[org].gulp(self.data[data])) 
                for org, data in tests ]

        print ' - Calculating ', len(todo), 'pairs...',
        sys.stdout.flush()

        results = p.map_async(gulp_call, todo)
        results.wait()
        results = results.get()
        print 'finished!'

        any_good = False
        killed = set()
        org_results = defaultdict(dict)
        for org, data, result in results:
            if not result:
                #print 'Killing', org
                self.kill(org)
                killed.add(org)
                continue
            org_results[org][data] = result

        if killed:
            print " - %s organisms didn't successfully evaluate" % len(killed)
        for k in killed:
            if k in org_results.keys():
                del org_results[k]

        self.fitness(org_results)
        self.generations[-1] = sorted(org_results.keys(), key=lambda k:
                self.organisms[k].fitness)

    def evaluate_generation(self, generation):
        self.bulk_evaluate(datas=self.fit_set,
                orgs=self.generations[generation])

    def evaluate_next(self):
        self.evaluate_generation(len(self.generations)-1)

    ### fitness function
    def fitness(self, results):
        org_ids = results.keys()
        data_ids = results.values()[0].keys()

        net_e_err = 0.0
        net_s_err = 0.0
        net_f_err = 0.0

        for o in org_ids:
            org = self.organisms[o]
            org.e_err = 0.0
            org.s_err = 0.0
            org.f_err = 0.0
            for d in data_ids:
                e_err = results[o][d]['energy'] - self.data[d].energy
                s_err = results[o][d]['stresses'] - self.data[d].stresses
                f_err = results[o][d]['forces'] - self.data[d].forces
                org.e_err += abs(e_err)
                org.s_err += np.average(abs(s_err))
                org.f_err += np.average(np.average(abs(f_err)))
            org.energy_err = org.e_err/len(data_ids)
            org.stress_err = org.s_err/len(data_ids)
            org.force_err = org.f_err/len(data_ids)
            net_e_err += org.energy_err
            net_s_err += org.stress_err
            net_f_err += org.force_err

        for o in org_ids:
            org = self.organisms[o]
            org.energy_fitness = org.energy_err/net_e_err
            org.stress_fitness = org.stress_err/net_s_err
            org.force_fitness = org.force_err/net_f_err
            org.fitness = sum([
                self.weights['energy']*org.energy_fitness,
                self.weights['stress']*org.stress_fitness,
                self.weights['force']*org.force_fitness])

    ### Selection operators

    def _tournament_select(self):
        orgs=[]
        for i in range(self.tourn_size):
            orgs.append(random.choice(self.generations[-1]))
        orgs = sorted(orgs, key=lambda k:
                self.organisms[k].fitness)
        return self.organisms[orgs[0]]

    def _roulette_select(self):
        total = sum([ self.organisms[k].fitness for k in self.generations[-1] ])
        score = random.random()
        orgs = sorted(self.generations[-1], key=lambda k:
                self.organisms[k].fitness)
        current = 0.0
        for k in orgs:
            current += self.organsims[k].fitness/total
            if current > score:
                return self.organisms[k]

    def select(self):
        return self._tournament_select()

    ### mating operators

    def _mix_mate(self, org1, org2):
        child = self.create_organism()
        for k in self.genome.cont_vars:
            p1_frac = random.random()
            child.genes[k] = p1_frac*org1.genes[k]+(1-p1_frac)*org2.genes[k]
        for k in self.genome.disc_vars:
            if random.random() < 0.5:
                child.genes[k] = org1[k]
            else:
                child.genes[k] = org2[k]
        return child

    def mate(self, org1, org2):
        child = self._mix_mate(org1, org2)
        while not child.valid:
            child = self._mix_mate(org1, org2)
        return child

    ### mutation operators

    def _mutate(self, org):
        if isinstance(org, int):
            org = self.organisms(org)
        for k in org.cont_vars:
            if random.random() < self.p_mutate:
                org.genes[k] = random.choice(
                        np.linspace(org.cont_vars[k][0], org.cont_vars[k][1], 1e3))
        for k in org.disc_vars:
            if random.random() < self.p_mutate:
                org.genes[k] = random.choice(org.disc_vars[k])
        return org

    def mutate(self, org):
        org = self._mutate(org)
        while not org.valid:
            self._mutate(org)
        return org

    ### Population management

    def create_organism(self):
        org = self.genome()
        org.id = len(self.organisms)
        self.organisms[org.id] = org
        return org

    def initialize_population(self):
        gen = []
        while len(gen) < self.first_generation_factor*self.pop_size:
            child = self.create_organism()
            gen.append(child.id)
        self.generations= [gen]

    def create_generation(self):
        gen = []
        for i in range(int(round(self.pop_size*(1-self.f_replace)))):
            gen.append(self.generations[-1][i])

        while len(gen) < self.pop_size:
            parent1 = self.select()
            parent2 = self.select()
            child = self.mate(parent1, parent2)
            child = self.mutate(child)
            gen.append(child.id)

        self.generations.append(gen)

    ### Optional optimization stuff

    def rescale(self):
        if len(self.organisms) < 1000:
            return
        for k, data in zip(self.genes, self.keys()):
            if k not in self.genome.cont_vars:
                continue


    def shift(self):
        if len(self.organisms) < 1000:
            return
        for k, data in zip(self.genes, self.keys()):
            if k in self.genome.cont_vars:
                self.genome.cont_vars[k] = prune_range(data)
            if k in self.genome.disc_vars:
                self.genome.disc_vars[k] = prune_values(data)

    def refine_genome(self):
        self.rescale()
        self.shift()

    ### Fit data selection / validation

    def _random_fit_data(self):
        self.fit_set = []
        while len(self.fit_set) < self.fit_size:
            self.fit_set.append(random.choice(self.data.keys()))

    def select_fit_data(self):
        return self._random_fit_data()

    ### data logging
    def initialize_outputs(self):
        self.keys = self.organisms[0].genes.keys()
        f = open('ipr.log','w')
        self.csvf = csv.DictWriter(f, self.keys)        

    def update_results(self):
        if not hasattr(self, 'csvf'):
            self.initialize_outputs()
        rows = [ self.organisms[k].genes for k in self.generations[-1] ]
        self.csvf.writerows(rows)

    def update_best(self):
        self.best.save('best.pot')

    def summarize(self):
        print ' - Best energy fitness:', self.best.energy_err
        print ' - Run time:', self.laps[-1]-self.laps[-2], ' seconds' 


    def output(self):
        self.update_best()
        self.update_results()
        self.summarize()

    ### Actual run command

    def __call__(self, generations=300):
        '''
        Fits a potential given the specified inputs. The manner in which the
        optimizer works can be modified by changing functions of the Optimizer
        instance you are running. Functions are pluggable with arbitrary
        replacements as long as the meet the requirements of the function's
        documentation.

        GA cycle:

        self.select_fit_data()
        self.initialize_population()
        self.evaluate_next()

        for gen in range(generations):
            self.create_generation()
            self.evaluate_next()
            self.output()
        '''
        print 'Generation 0:\n================'
        print ' - Initializing'

        self.laps = [time.time()]
        self.select_fit_data()
        self.initialize_population()
        self.evaluate_next()
        self.laps.append(time.time())

        for gen in range(generations):
            print 'Generation %s:\n================' % len(self.generations)
            self.create_generation()
            self.evaluate_next()
            self.laps.append(time.time())
            self.output()
            #self.refine_genome()
