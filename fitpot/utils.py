from data import *
import numpy as np
import os
import random

# Function to standardize data (center about 0, divide by stdev to eliminate units)
def normalize(data):
    data = np.array(data)
    return data - np.average(data)/np.std(data)

# Function to normalize data (all values in the range [0,1])
def scale(data):
    data = np.array(data)
    minval = min(data)
    maxval = max(data)
    return (data-minval)/(maxval-minval)

def unitize(data):
    data = np.array(data)
    return data/sum(data)

def read_frcout(frcout):
    energy = 0
    forces = []
    stresses = []
    do_forces = False
    do_stresses = False
    for line in frcout:
        if 'energy' in line:
            try:
                energy = float(line.split()[1])
            except ValueError:
                raise FortranError
        if 'gradients cartesian' in line:
            do_forces = True
            continue
        if 'strain' in line:
            do_forces = False
            do_stresses = True
            continue
        if do_forces:
            force = []
            if len(line.split()) != 4:
                raise FortranError
            for f in line.split()[1:4]:
                try:
                    force.append(float(f))
                except ValueError:
                    raise FortranError
            forces.append(force)
        if do_stresses:
            if len(line.split()) != 3:
                raise FortranError
            for s in line.split():
                try:
                    stresses.append(float(s))
                except ValueError:
                    raise FortranError
    return {'energy':energy,
            'forces':np.array(forces),
            'stresses':np.array(stresses)}

def ensure_length(line):
    if len(line) > 80:
        llist = line.split()
        mid = len(llist)/2
        l1 = ' '.join(llist[:mid+1])
        l2 = ' '.join(llist[mid+1:])
        return '%s &\n %s' % (ensure_length(l1), ensure_length(l2))
    return line

