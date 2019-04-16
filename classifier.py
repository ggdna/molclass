import numpy as np
import pandas as pd
from sklearn import svm

lscale = 1

def acquisition(datafile, p2gfile):
    '''Takes files in form of classdata.csv and probe_to_gene.csv and returns a
       dataframe that composes the two'''
    data = pd.read_csv(datafile)
    genes = pd.read_csv(p2gfile)
    raw = data.merge(genes, how='left', left_on='ID', right_on='ID')
    return raw

def dftomatrix(df, m):
    '''Takes training or testing set as a dataframe and number of rows and
       converts to 1D matrix of values (no labels)'''
    ndf = df[:m]
    dfv = ndf.iloc[:,2:].values.tolist()
    dfv = list(zip(*dfv))
    return dfv

def categorize(raw, param1, param2, n):
    '''Takes raw data set delimeters for classes (param1 and param2) and n for
       which the testing set is 1/n of total size and number of rows of data to
       consider, and outputs curated training and testing value and label matrices'''
    c1train = []
    c1test = []
    c2train = []
    c2test = []
    c1c = 0
    c2c = 0
    param1 = param1.split()
    param2 = param2.split()
    for col in raw.columns:
        if any(p in col for p in param1):
            c1c += 1
            if (c1c%n)==0:
                c1test.append(col)
            else:
                c1train.append(col)
        elif any(p in col for p in param2):
            c2c += 1
            if (c2c%n)==0:
                c2test.append(col)
            else:
                c2train.append(col)
    dftrain = pd.concat([raw.ID, raw.Gene, raw[c1train], \
                         raw[c2train]], axis=1, join='inner')
    dftest = pd.concat([raw.ID, raw.Gene, raw[c1test], \
                        raw[c2test]], axis=1, join='inner')
    if len(c1test) == 0:
        dftest = dftrain
        c1test = c1train
        c2test = c2train
    trainvals = dftomatrix(dftrain, len(raw))
    testvals = dftomatrix(dftest, len(raw))
    if lscale:
        trainvals = [tuple(2**i for i in tval) for tval in trainvals]
        testvals = [tuple(2**i for i in tval) for tval in testvals]
    trainlabels = [1]*len(c1train) + [0]*len(c2train)
    testlabels = [1]*len(c1test) + [0]*len(c2test)
    return trainvals, testvals, trainlabels, testlabels

def train(trainvals, trainlabels, C, class_weight):
    '''For given 1D training matrices and classweight dict and lambda values, using
       sklearn, an SVM will be trained, and weight matrix will be returned'''
    clf = svm.LinearSVC(loss='squared_hinge', penalty='l1', dual=False,
                        C=C, class_weight=class_weight)
    clf.fit(trainvals, trainlabels)
    weights = clf.coef_.flatten().tolist()
    thresh = clf.intercept_[0]
    if thresh != 0:
        print("Threshold not 0.")
        print(thresh)
    return weights

def normalize(weights, n):
    '''Normalizes weight matrix to desired value and rounds to nearest whole'''
    mags = [abs(w) for w in weights]
    mw = max(mags)
    if mw == 0:
        norm_factor = 0
    else:
        norm_factor = float(n)/mw
    norm_weights = [np.around(w*norm_factor) for w in weights]
    return norm_weights

def subsetsums(l):
    '''Returns the sums of all subsets of a list'''
    def combinations(seq):
        yield seq
        for i in range(len(seq)):
            for combination in combinations(seq[:i] + seq[i+1:]):
                yield combination
    return map(sum, combinations(l))

def realize(weights, m):
    '''Returns implementable approximation of desired weights from simulated probe
       behavior.'''
    rweights = {}
    nzweights = {}
    probes = {}
    for i in range(len(weights)):
        if weights[i] != 0:
            probes[i] = [np.sign(weights[i])*10**n for n in list(np.random.normal(0,0.5,m))]
            nzweights[weights[i]] = i
    mw = sum(probes[nzweights[sorted(nzweights.iterkeys())[0]]])
    normf = mw/(sorted(nzweights.iterkeys())[0])
    rweights[nzweights[sorted(nzweights.iterkeys())[0]]] = mw
    for w in sorted(nzweights.iterkeys())[1:]:
        rweights[nzweights[w]] = subsetsums(probes[nzweights[w]])
        rweights[nzweights[w]] = min(rweights[nzweights[w]], key=lambda x:abs(x-(w*normf)))
    return [rweights[i] if i in rweights.keys() else weights[i] for i in range(len(weights))]

def targetgenes(weights, raw):
    '''For given weight matrix and raw dataframe, the targeted genes (the genes for
       which weights are nonzero) are determined and output with their respective
       weights in a dictionary'''
    output = {}
    for i in range(len(weights)):
        if weights[i] != 0:
            output[raw.Gene[i]] = weights[i]
    return output

def test(weights, testvals, testlabels):
    '''For given weight matrix and testing set, the accuracy of the SVM is
       evaluated as a percentage, and returned with the size of the testing set'''
    test_res = [np.dot(tval, weights) for tval in testvals]
    correct = 0
    testsize = len(test_res)
    for i in range(testsize):
        if test_res[i] > 0:
            if testlabels[i] == 1:
                correct += 1
        else:
            if testlabels[i] == 0:
                correct += 1
    accuracy = float(correct)/testsize
    return accuracy, testsize

def mutate(mu, sig, wmatrix):
    '''For given weight matrix and mean and standard deviation (of exponent by which
        weight is modified), a mutated weight matrix is returned'''
    nwmatrix = []
    for i in range(len(wmatrix)):
        if wmatrix[i] != 0:
            nwmatrix.append(10**(np.random.normal(mu, sig, 1)[0])*wmatrix[i])
        else:
            nwmatrix.append(0)
    return nwmatrix

def lambdaopt(trainvals, trainlabels, testvals, testlabels, m0, m1, n, \
              maxfeats, class_weight={0:1, 1:1}):
    '''For testing set, training set, normalization parameter, and min/max/steps
        (m0,m1,n) for lambda value search space, and maximum desired number of
        features lambda value and corresponding weight matrix are returned,
        with accuracy'''
    lams = np.linspace(m0, m1, n)
    topweights = []
    topacc = 0
    topfeat = len(trainvals[0])
    toplam = lams[0]
    for lam in lams:
        weights = train(trainvals, trainlabels, lam, class_weight)
        print("Trained!")
        nonzero = [w for w in weights if w != 0]
        if len(nonzero)<=maxfeats:
            accuracy = test(weights, testvals, testlabels)[0]
            if accuracy>topacc or (accuracy == topacc and len(nonzero)<topfeat):
                topacc = accuracy
                topweights = weights
                topfeat = len(nonzero)
                toplam = lam
        print(lam)
    return topweights, toplam, topacc, topfeat
