import numpy as np
from numpy.linalg import inv
import os, os.path

### IPR specific stuff ###
class FortranError(Exception):
    '''When fortran doesn't keep spaces between numbers so everything goes to
    shit because fortran is such a decrepit POS'''

class GulpError(Exception):
    '''gulp runtime error'''

class VaspError(Exception):
    '''gulp runtime error'''

class DataError(Exception):
    '''Incorrectly assembled snapshot'''

### IPR specific stuff ###

class Data:
    def __init__(self):
        '''
        Data objects must have an input, and output attributes.

        inputs:
            cell        = [[1, 0, 0], [0, 1, 0], [0, 0, 1]] 
            atoms       = [ ('A', [0,0,0], 'B', [0.5,0.5,0.5]) ]
        outputs:
            energy      = float()
            forces      = [[0, 0, 0], [0, 0, 0]]
            stresses    = [0, 0, 0, 0, 0, 0]
        '''
        self.cell = [[]]
        self.atoms = []
        self.energy = 0.0
        self.stresses = [0,0,0,0,0,0]
        self.forces = [[]]

    @staticmethod
    def read_path(search_path):
        data = []
        for (path, dirs, files) in os.walk(search_path):
            for file in files:
                if 'OUTCAR' in file:
                    data += Data.read_outcar(path+'/'+file)
            if set(['energy', 'stresses', 'forces', 'POSCAR']) < set(files):
                data.append(Data.read_snapshot(path))
        return data

    @property
    def inputs(self):
        return {'cell':self.cell, 
                'atoms':self.atoms}

    @property
    def outputs(self):
        return {'energy': self.energy,
                'forces': self.forces,
                'stresses': self.stresses}
        
    @staticmethod 
    def read_snapshot(snapdir):
        ### test that path is a snapshot ###
        if not ( os.path.exists(snapdir+'/POSCAR') and
                os.path.exists(snapdir+'/energy') and
                os.path.exists(snapdir+'/forces') and
                os.path.exists(snapdir+'/stresses')):
            raise DataError
        
        snapshot = Data()
        snapshot.energy = float(open(snapdir+'/energy').read())
        snapshot.stresses = [ float(f) for f in 
                open(snapdir+'/stresses').read().split() ]
        poscar = open(snapdir+'/POSCAR').readlines()
        snapshot.cell = [ [ float(f) for f in vec.split() ] 
                for vec in poscar[2:6] ]

        ## test vasp4/5, direct/cart
        atom_counts = []
        v5line = poscar[6].lower()
        if v5line[0].isalpha():
            ### has a species line -> vasp5
            atom_types = v5line.split()
            atom_types = [ d.split('_')[0] for d in poscar ]
            atom_counts = [ int(f) for f in poscar[7].split() ]
            dline = poscar[8].lower()
            atoms = poscar[9:sum(atom_counts)+1]
        else:
            ### doesn't -> vasp4
            if not os.path.exists(snapdir+'/atom_types'):
                raise DataError
            atom_types = open(snapdir+'/atom_types').read()
            atom_types = atom_types.strip().split()
            atom_counts = [ int(f) for f in dline ]
            v5line = poscar[7].lower()
            atoms = poscar[8:sum(atom_counts)+1]

        atom_array = []
        for type, count in zip(atom_types, atom_counts):
            atom_array += [type]*count

        if v5line.startswith('d'):
            cart = False
        else:
            cart = True
            inv_cell = inv(snapshot.cell)

        for elt, atom in zip(atom_array, atoms):
            coord = [ float(f) for f in atom.split() ]
            if cart:
                coord = dot(inv_cell, coord)
            snapshot.coords.append((elt, list(coord)))

        forces = open(snapdir+'/forces').readlines()
        if not len(forces) == len(snapshot.coords):
            raise DataError
        snapshot.forces = [ [ float(f) for f in line ]  for line in forces ]
        
        return snapshot

    @staticmethod
    def read_outcar(outcar):
        data = open(outcar).readlines()
        snapshots = []
        atom_types = []
        atom_counts = []
        atom_array = []

        snapshot = Data()

        for n,line in enumerate(data):
            if 'POTCAR:' in line:
                temp = line.split()[2]
                for c in ['.','_','1']:
                    if c in temp:
                        temp = temp[0:temp.find(c)]
                atom_types += [temp]
            if 'ions per type' in line:
                atom_types = atom_types[:len(atom_types)/2]
                atom_counts = [ int(f) for f in line.split()[4:] ]
                for type, count in zip(atom_types, atom_counts):
                    atom_array += [type]*count
            if 'direct lattice vectors' in line:
                cell = []
                for i in range(3):
                    temp = data[n+1+i].split()
                    cell += [[float(temp[0]), float(temp[1]), float(temp[2])]]
                inv_cell = inv(cell)
                snapshot.cell = cell
            if 'FREE ENERGIE OF THE ION-ELECTRON SYSTEM' in line:
                snapshot.energy = float(data[n+4].split()[6])
            if 'STRESS in cart' in line:
                for iline in range(20):
                    nline = data[n+iline]
                    if 'Total' in nline:
                        snapshot.stresses = [ float(f) for f in 
                                nline.split()[1:]]
                        break
            if 'POSITION          ' in line:
                forces = []
                atoms = []
                for iatom, elt in enumerate(atom_array):
                    temp = data[n+2+iatom].split()
                    forces += [[float(f) for f in temp[3:6]]]
                    atoms += [(elt, 
                        np.dot(inv_cell, [float(f) for f in temp[0:3]])
                        )]
                snapshot.forces = forces
                snapshot.coords = atoms
                snapshots.append(snapshot)
                snapshot = Data()
        return snapshots

