import pandas as pd
import secstruc as ss
import melting

r = 'gatatgcttcaggtcaggga'.upper()
lp = 'ccccgccaatatgagaaagtgtgtagcgcgatcagatgtt'.upper()
ln = 'gggattcatcggtagtagcagtgtagcgcgatcagatgtt'.upper()

def acquisition(sf):
    con = pd.read_csv(sf)
    return list(con["Gene"]), [s.lower() for s in con["Sequence"]]

def gc(seq):
    return float((seq.count('G')+seq.count('C')))/len(seq)

def sec(s1, s2):
    '''bigger is better'''
    pplp = ss.pairs([r+s1, s2+lp], material="dna")
    ppln = ss.pairs([r+s1, s2+ln], material="dna")
    pplp = [p[2] for p in pplp if p[1]==(len(r)+len(lp)+len(s1)+len(s2)+1)]
    ppln = [p[2] for p in ppln if p[1]==(len(r)+len(lp)+len(s1)+len(s2)+1)]
    return min(pplp)+min(ppln)

def score(seq):
    s1, s2 = seq[int(len(seq)/2):], seq[:int(len(seq)/2)]
    if abs(melting.temp(s1)-72.5)<7.5 and abs(melting.temp(s2)-72.5)<7.5 and abs(gc(s1)-0.5)<0.1 and abs(gc(s2)-0.5) and (s2[0]=="A" or s2[0]=="T"):
        ssc = sec(s1, s2)
        return ssc
    return 0

def pick(tr):
    scores = {}
    for i in range(len(tr)-40):
        seq = tr[i:(i+40)].upper()
        scores[score(seq)] = seq
    return scores[max(scores.keys())]

def probes(sf):
    probes = []
    genes, seqs = acquisition(sf)
    for i in range(len(genes)):
        s = pick(seqs[i])
        s1, s2 = s[int(len(s)/2):], s[:int(len(s)/2)]
        probes.append([genes[i], r+s1, s2+lp, s2+ln])
    return probes
