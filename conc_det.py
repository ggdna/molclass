import numpy as np
import pandas as pd

def acquisition(pbf):
    data = pd.read_csv(pbf)
    genes = list(data["Gene"])
    weights = list(data["Weight"])
    bias = [data["Reads"][i]/data["Expression"][i] for i in range(len(weights))]
    return genes, weights, bias

def mixingconc(weights, bias):
    '''Given a list of weights and respective biases, a list of mixing
       concentration ratios is determined, given as a list of [posconc, negconc]'''
    conc = [None]*len(weights)
    normf = None
    for j in range(len(weights)):
        if weights[j] != 0:
            normf = abs(bias[j]/weights[j])
            break
    i = 0
    while i != len(weights):
        if weights[i] != 0:
            k = (0.5*normf*weights[i])/bias[i]
            conc[i] = [0.5+k,0.5-k]
            if np.around(abs(k), 5) > 0.5:
                normf = abs(bias[i]/weights[i])
                i = 0
            else:
                i += 1
        else:
            conc[i] = [0,0]
            i += 1
    return conc

def conctoweights(conc, bias, prec):
    '''Given concentrations of positive/negative weights, determine effective weights'''
    weights = [0]*len(conc)
    for i in range(len(conc)):
        weights[i] = 2*(np.around(conc[i][0], decimals=prec)-0.5)*bias[i]
    return weights
