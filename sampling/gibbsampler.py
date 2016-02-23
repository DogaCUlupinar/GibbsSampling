'''
Created on Feb 20, 2016

@author: dulupinar
'''
from collections import defaultdict,Counter
import operator
import random
import numpy as np
import logging
logging.basicConfig(level=logging.CRITICAL,format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')


ALPHA = ["A","C","T","G"]
LEN_ALPHA = len(ALPHA)

def generateSequences(genome_length,kmer_length,num_sequneces,mismatches,rand_seed=1):
    #check that genome_length is greater than kmer_length
    random.seed(rand_seed)
    kmer_no_mismatch = [ALPHA[random.randint(0,LEN_ALPHA-1)] for i in range(kmer_length)]
    sequences = []
    for i in range(num_sequneces):
                
        genome = [ALPHA[random.randint(0,LEN_ALPHA-1)] for i in range(genome_length-kmer_length)]
        kmer = list(kmer_no_mismatch)
        #generate the mismatches
        for _ in range(random.randint(0,mismatches)):
            # up to mismatches
            rand_pos_kmer = random.randint(0,kmer_length-1) #random position in kmer
            kmer[rand_pos_kmer] = ALPHA[(ALPHA.index(kmer[rand_pos_kmer])+1)%4]
        
        kmer = ''.join(kmer)
        logging.debug("inserted kmer {0} with {1} mismatches".format(kmer,mismatches))
        #put it in the genome
        rand_pos_genome = random.randint(0,genome_length-kmer_length -1)
        genome.insert(rand_pos_genome, kmer)
        sequences.append("".join(genome))
        
    return sequences,''.join(kmer_no_mismatch)

class SequenceKmerList():
    
    def __init__(self,sequence_list,kmer_length):
        self.kmer_length = kmer_length
        self.sequences_kmer = [[sequence,""] for sequence in sequence_list]
        
    def getAllButOneKmers(self):
        #doest return one kmer that is the excluded index
        num_sequences = len(self.sequences_kmer)
        exclude_index = random.randint(0,num_sequences -1)
        return [sequence_kmer[1] for index,sequence_kmer in enumerate(self.sequences_kmer) if index != exclude_index],exclude_index
    
    def __getitem__(self,index):
        return self.sequences_kmer[index]
        
    def __len__(self):
        return len(self.sequences_kmer)
    
    def getAllKmers(self):
        return [sequence_kmer[1] for sequence_kmer in self.sequences_kmer]

def calculateScore(kmers):
    #calculates the score given a list of kmers
    #calucluate a column for each
    
    #create the columns
    columns = [ Counter() for _ in kmers[0]]
    for kmer in kmers:
        for index,base in enumerate(kmer):
            columns[index][base]+=1
    
    score = 0
    for column in columns:
        most_common_in_column = column.most_common(1)[0][0]
        for base in ALPHA:
            if base != most_common_in_column:
                score+=column[base]
    
    consensus = "".join([col.most_common(1)[0][0] for col in columns])
    return score,consensus
    

def gibbs(sequences,kmer_length,rand_seed=1):
    random.seed(rand_seed)
    sequence_kmers = SequenceKmerList(sequences,kmer_length)  
    
    logging.debug("Selecting initial random kmers for each sequence")
    for i in range(len(sequences)):        
        sequence = sequence_kmers[i][0]
        start_pos = random.randint(0,len(sequence) - kmer_length)
        end_pos = start_pos+ kmer_length
        kmer = sequence[start_pos: end_pos]
        sequence_kmers[i][1] = kmer
        
    score = float('Inf')
    tries = 0
    max_tries = 10000
    max_changless = 100 # number of tries we want to do without score changing
    score_change = 0 #number of times where the score hasn't changed
    min_score = float('Inf')
    while score > 1 and tries < max_tries and score_change<max_changless:
        tries +=1
        if tries % (max_tries/15) == 0: 
            logging.debug("Score is " + str(score) + " on try " + str(tries))
        
        gibbsSampler(sequence_kmers)
        score,consensus = calculateScore(sequence_kmers.getAllKmers())
        
        if score < min_score:
            score_change = 0
            min_score = score
            min_consensus = consensus
            logging.debug("found min score {0} and consensus {1} on try {2}".format(min_score,min_consensus,tries))
        else:
            score_change +=1    

    logging.warn("The max tries is {0} and changless is {1}".format(tries,score_change))
    return min_consensus,min_score
        
def gibbsSampler(sequence_kmer_list):
    
    kmer_list,exclude_index = sequence_kmer_list.getAllButOneKmers()
    profile = generateProfile(kmer_list) 
    excluded_sequence = sequence_kmer_list[exclude_index][0]   
    #select kmer
    kmer_selection = dict()
    sum_prob = 0
    for i in range(len(excluded_sequence) - sequence_kmer_list.kmer_length +1):
        possible_kmer = excluded_sequence[i:i+sequence_kmer_list.kmer_length]
        #generate probability for kmer
        probability = 1
        for index,base in enumerate(possible_kmer):
            probability = probability*profile[ALPHA.index(base),index]
        kmer_selection[possible_kmer] = probability
        sum_prob+=probability
    #roll the dice need to roll
    prev = 0
    for key in kmer_selection:
        kmer_selection[key] =  kmer_selection[key]/sum_prob + prev
        prev = kmer_selection[key]
        
    dice_roll = random.uniform(0,1)
    for value in sorted(kmer_selection.items(), key=operator.itemgetter(1)):
        if value[1] > dice_roll:
            break
    selected_kmer = value[0]
    
    sequence_kmer_list[exclude_index][1] = selected_kmer
    
    return sequence_kmer_list

def generateProfile(kmers):  
    #select the kme
    
    profile = np.zeros((LEN_ALPHA,len(kmers[0]))) #profile is ACTG by count at position
    for kmer in kmers:
        for index,base in enumerate(kmer):
            profile[ALPHA.index(base),index]+=1
    
    #add laplace smoothing
    for cell in np.nditer(profile,op_flags=['readwrite']):
        cell[...] = cell +1
        
    return profile
 
def runRegularGibbs():           
    seqs,kmer = generateSequences(25, 5, 10, 1)
    print "THE ORIGINAL SEQUENCES"
    print str("\n".join(seqs))
    print "THE ORIGINAL KMER"
    print kmer
 
    found_kmer = defaultdict(list)
    for i in range(100):
        km = gibbs(seqs, len(kmer),rand_seed =None )
        found_kmer[km[1]].append(km[0])
        
    ks = min(found_kmer.keys())
    print ks,found_kmer[ks]
    print calculateScore(found_kmer[ks])
       
if __name__ == "__main__":
    
    seqs,kmer = generateSequences(35, 15, 10, 3)
    print "THE ORIGINAL SEQUENCES"
    print str("\n".join(seqs))
    print "THE ORIGINAL KMER"
    print kmer
 
    found_kmer = defaultdict(list)
    for i in range(100):
        km = gibbs(seqs, len(kmer),rand_seed =None )
        found_kmer[km[1]].append(km[0])
        
    ks = min(found_kmer.keys())
    print ks,found_kmer[ks]
    print calculateScore(found_kmer[ks])