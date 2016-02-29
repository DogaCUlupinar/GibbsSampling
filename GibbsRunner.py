from sampling.gibbsampler import *
from collections import defaultdict,Counter
from operator import itemgetter
import random
import numpy as np
import logging
import math
import matplotlib.pyplot as plt
import logging

CONFIG_LABEL  = "LABEL"
CONFIG_PRIOR  = "PRIOR"
CONFIG_GREEDY = "GREEDY"
CONFIG_COLOR  = "COLOR"
configurations = [{CONFIG_LABEL:'Standard Gibbs', CONFIG_GREEDY:False,CONFIG_PRIOR:False, CONFIG_COLOR:'b' },\
                   {CONFIG_LABEL:'Prior Knowledge', CONFIG_GREEDY:False,CONFIG_PRIOR:True,CONFIG_COLOR:'r' },\
                   {CONFIG_LABEL:'Greedy', CONFIG_GREEDY:True,CONFIG_PRIOR:False,CONFIG_COLOR:'g'},\
                   {CONFIG_LABEL:'Prior Knowledge and Greedy', CONFIG_GREEDY:True,CONFIG_PRIOR:False, CONFIG_COLOR:'k' }]
    
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
    return ks,found_kmer[ks]

def makeNumIterations(seqs,kmer_len,ax,configurations,passes,original_kmer):           

    max_iterations = 1000
    found_kmer = defaultdict(list)
    count= Counter()
    
    sequence_kmers = selectInitalKmers(seqs, kmer_len)
    
    total_iterations = []
    for pass_c in range(passes):
        if pass_c % max(1,(passes/10)) == 0: logging.warn("{0:.0f}% is done".format(float(pass_c)/passes*100))
        
        for i in range(max_iterations):  
            consensus_string,score,data,profile = gibbs(sequence_kmers,determinstic=configurations[CONFIG_GREEDY] )
            prior_profile = profile if configurations[CONFIG_PRIOR] else None
            sequence_kmers = selectInitalKmers(seqs, kmer_len,prior = prior_profile)
            count[score]+=1
            found_kmer[score].append(consensus_string)
            if consensus_string == original_kmer:
                break
        total_iterations.append(i)
    
    lables,values = zip(*sorted(count.items(),key=itemgetter(0)))
    ax.plot(lables,values,'-o',label=configurations['LABEL'])
    
    print configurations[CONFIG_LABEL]
    print "Average Iteration till motif: ",sum(total_iterations)/float(passes),np.std(total_iterations)
    ks = min(found_kmer.keys())
    return ks,found_kmer[ks]

def makeSTDevGraph(seqs,kmer_len,ax,configurations,iterations):           
    
    found_kmer = defaultdict(list)
    count= Counter()
    lines = []
    iter_to_conv = [] #iterations required to converge
    
    sequence_kmers = selectInitalKmers(seqs, kmer_len)
    
    for i in range(iterations):
        if i % max(1,(iterations/10)) == 0: print "{0:.0f}% is done".format(float(i)/iterations*100)
        consensus_string,score,data,profile = gibbs(sequence_kmers,determinstic=configurations[CONFIG_GREEDY] )
        prior_profile = profile if configurations[CONFIG_PRIOR] else None
        sequence_kmers = selectInitalKmers(seqs, kmer_len,prior = prior_profile)
        count[score]+=1
        found_kmer[score].append(consensus_string)
        lines.append(data)
        iter_to_conv.append(data[0][-1])
    
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
    
    ax.plot(mean_curve,'-o',label=configurations[CONFIG_LABEL],color=configurations[CONFIG_COLOR])
    ax.errorbar(np.arange(max_iteration),mean_curve,yerr=stdev_curve,color=configurations[CONFIG_COLOR], linestyle="None")
    
    ks = min(found_kmer.keys())
    print configurations[CONFIG_LABEL]
    print "Average iteration is {0}, with sdt {1}".format(str(sum(iter_to_conv)/float(len(iter_to_conv))),str(np.std(np.asarray(iter_to_conv))))
    return ks,found_kmer[ks]

def makeSTDGraphs(mismatch,iter):
    fig,ax=plt.subplots()
    mismatch = 0
    seqs,kmer = generateSequences(250, 15, 10,mismatch)
    print "THE ORIGINAL SEQUENCES"
    print str("\n".join(seqs))
    print "THE ORIGINAL KMER"
    print kmer
    print "THE NUMBER OF MISMATCHES ",mismatch
    
    for config in configurations:
        print makeSTDevGraph(seqs,len(kmer),ax,config,iter)
    
    
    ax.legend(loc='upper right')
    plt.xlabel('Iterations')
    plt.ylabel('Motif Entropy')
    plt.title("Motifs with {0} mismatch(es)".format(str(mismatch)))
    plt.show()

def makeIterGraphs(mismatch,iter):
    fig,ax=plt.subplots()
    seqs,kmer = generateSequences(250, 15, 10, mismatch)
    print "THE ORIGINAL SEQUENCES"
    print str("\n".join(seqs))
    print "THE ORIGINAL KMER"
    print kmer
    print "THE NUMBER OF MISMATCHES ",mismatch
    

    for config in configurations:
        print makeNumIterations(seqs,len(kmer),ax,config,iter,kmer)
    
    ax.legend(loc='upper right')
    plt.xlabel('Motif Entropy')
    plt.ylabel('Count')
    plt.title("Motifs with {0} mismatch(es)".format(str(mismatch)))
    plt.show()
    
if __name__ == "__main__":
    makeIterGraphs(3,1000)