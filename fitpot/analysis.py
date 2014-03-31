import numpy as np
import matplotlib.pylab as plt
from collections import defaultdict

def gene_plot(pfile):
    data = get_data(pfile)
    for gene in data:
        res = plt.hist(gene)
        plt.show()

def get_data(pfile):
    f = open(pfile)
    data = [ [ float(ff) for ff in d.split(',') ] for d in f.readlines() ]
    data = np.array(data)
    return data.T

#### intended for use within optimizer

def data_at_extremum(data):
    hist, edges = np.histogram(data)
    tot = sum(hist)
    nhist = [ float(h)/tot for h in hist ]
    if nhist[-1] > 0.5 or nhist[0] > 0.5:
        return True
    else:
        return False

def prune_range(data):
    hist, edges = np.histogram(data)
    tot = sum(hist)
    nhist = [ float(h)/tot for h in hist ]
    max_val = max(nhist)
    min, max = edges[0], edges[-1]
    while nhist[0] < 1e-2:
        nhist.pop(0)
        edges.pop(0)

    while nhist[-1] < 1e-2:
        nhist.pop(-1)
        edges.pop(-1)
    return [edges[0], edges[-1]]

def prune_values(data, threshold=0.05):
    counts = defaultdict(int)
    tot = len(data)
    for d in data:
        counts[d] += 1
    scounts = dict( (k, v/float(tot)) for k, v in counts )
    return [ k for k, v in scounts.values() if v > 0.5 ]

def shift_range(data, threshold=0.25):
    hist, edges = np.histogram(data)
    step = edges[0]-edges[1]
    nhist = [ float(h)/tot for h in hist ]
    bot = edges[0]
    top = edges[-1]

    if nhist[0] > threshold:
        bot -= step
    if nhist[-1] > threshold:
        top += step

    return [top, bot]
