from sampling.gibbsampler import *
from collections import defaultdict,Counter
from operator import itemgetter
import random
import numpy as np
import logging
import math
import matplotlib.pyplot as plt

CONFIG_LABEL  = "LABEL"
CONFIG_PRIOR  = "PRIOR"
CONFIG_GREEDY = "GREEDY"
configurations = [{CONFIG_LABEL:'Standard Gibbs', CONFIG_GREEDY:False,CONFIG_PRIOR:False },\
                   {CONFIG_LABEL:'Prior Knowledge', CONFIG_GREEDY:False,CONFIG_PRIOR:True },\
                   {CONFIG_LABEL:'Greedy', CONFIG_GREEDY:True,CONFIG_PRIOR:False },\
                   {CONFIG_LABEL:'Prior Knowledge and Greedy', CONFIG_GREEDY:True,CONFIG_PRIOR:False }]
    
def makeMotifScoreCountGraph(seqs,kmer_len,ax,configurations,iterations):           

 
    found_kmer = defaultdict(list)
    count= Counter()
    
    sequence_kmers = selectInitalKmers(seqs, kmer_len)
    
    for i in range(iterations):
        if i % max(1,(iterations/10)) == 0: print "{0:.0f}% is done".format(float(i)/iterations*100)
        consensus_string,score,data,profile = gibbs(sequence_kmers,determinstic=configurations[CONFIG_GREEDY] )
        prior_profile = profile if configurations[CONFIG_PRIOR] else None
        sequence_kmers = selectInitalKmers(seqs, kmer_len,prior = prior_profile)
        count[score]+=1
        found_kmer[score].append(consensus_string)
        
    lables,values = zip(*sorted(count.items(),key=itemgetter(0)))
    ax.plot(lables,values,'-o',label=configurations['LABEL'])
    
    ks = min(found_kmer.keys())
    return calculateScore(found_kmer[ks])
    
def makeSTDevGraph(seqs,kmer_len,ax,configurations,iterations):           
    
    found_kmer = defaultdict(list)
    count= Counter()
    lines = []
    
    sequence_kmers = selectInitalKmers(seqs, kmer_len)
    
    for i in range(iterations):
        if i % max(1,(iterations/10)) == 0: print "{0:.0f}% is done".format(float(i)/iterations*100)
        consensus_string,score,data,profile = gibbs(sequence_kmers,determinstic=configurations[CONFIG_GREEDY] )
        prior_profile = profile if configurations[CONFIG_PRIOR] else None
        sequence_kmers = selectInitalKmers(seqs, kmer_len,prior = prior_profile)
        count[score]+=1
        found_kmer[score].append(consensus_string)
        lines.append(data)
    
    """   
    for line in lines:
        ax.plot(*line,color='b')
    """
    
    max_iteration = reduce(lambda x,y : max(x,len(y[0])) if type(x) == int else max(len(x[0]),len(y[0])), lines) 
    mean_curve  = np.zeros(max_iteration)
    stdev_curve = np.zeros(max_iteration)
    
    #calculate mean and stdev     
    for i in range(max_iteration):
        col = []
        for line in lines:
            try:
                col.append(line[1][i])
            except IndexError:
                pass
        stdev_curve[i] = np.std(col)
        mean_curve[i] = sum(col)/float(len(col))
    
    ax.plot(mean_curve,'-o',label=configurations['LABEL'])
    ax.errorbar(np.arange(max_iteration),mean_curve,yerr=stdev_curve, linestyle="None")
    
    ks = min(found_kmer.keys())
    print ks,found_kmer[ks]
    print calculateScore(found_kmer[ks])
    return calculateScore(found_kmer[ks])

def makeGraphs():
    fig,ax=plt.subplots()
    seqs,kmer = generateSequences(250, 15, 10, 3)
    print "THE ORIGINAL SEQUENCES"
    print str("\n".join(seqs))
    print "THE ORIGINAL KMER"
    print kmer
    
    iter = 1000
    """
    for config in configurations:
        print makeMotifScoreCountGraph(seqs,len(kmer),ax,config,iter)
    """ 
    
    for config in configurations:
        print makeSTDevGraph(seqs,len(kmer),ax,config,iter)
    
    
    ax.legend(loc='upper right')
    plt.xlabel('Iterations')
    plt.ylabel('Motif Entropy')
    plt.show()     
if __name__ == "__main__":
    