'''
Created on Feb 23, 2016

@author: dulupinar
'''
import unittest
from gibbsampler import *
import matplotlib.pyplot as plt
import numpy as np

class TestSampler(unittest.TestCase):
    def test_no_mismatch(self):
        seqs,kmer = generateSequences(35, 15, 10, 0)
        self.assertEqual(gibbs(seqs, len(kmer) ),('GCCCTGAAGCGACAT', 13))
        
    def test_with_mismatch(self):
        seqs,kmer = generateSequences(35, 15, 10, 3)
        self.assertEqual(gibbs(seqs, len(kmer) ),('GCCCTGAAGCGACCA', 21))
        
    def test_slope(self):
        scores = [88, 87, 83, 81, 76, 68, 59, 48, 40, 32, 24, 13]
        print np.diff(scores)
        print scores
    
    def test_entropy(self):
        kmers = ["GTCG","GCTG","GCCT","CCCG","GGCG"]
        print calculateEntropy(generateProfile(kmers))
    
if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSampler)
    unittest.TextTestRunner().run(suite)