#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  3 12:36:26 2018

@author: grant
"""

import numpy as np
import time
from phase_genotype import phaseGenotype
import sys


def load_genome(file_name):
    gen = []
    f = open(file_name, 'r')
    for line in f:
        n = list()
        for i in line.strip():
            if i != " ":
                n.append(i)
        gen.append(n)
    f.close()
    #need to rearrange
    gen = np.array(gen, dtype=np.int8)
    gen = gen.transpose()
    return gen


def load_genome_faster(file_path):
    fileInput = open(file_path, 'r')
    firstLine = fileInput.readline()  # Manually get first line to detect size
    firstLine = map(int, firstLine.strip('\n').split(' '))  # Convert to list format
    genotypes = [ [] for individual in range(0, len(firstLine)) ]
    for individual in range(0, len(firstLine)):
        genotypes[individual].append(firstLine[individual])
    for line in fileInput:
        line = map(int, line.strip('\n').split(' '))
        for individual in range(0, len(line)):
            genotypes[individual].append(line[individual])
    fileInput.close()
    return genotypes



if __name__ == "__main__":
    
    filePath = "./CM124 Programming Assignment Guidelines/example_data_1.txt"
    #filePath = "TestData.txt"
    print "---------- Loading genotypes from file into memory ----------"
    start = time.time()
    genotypes = load_genome_faster(filePath)
    end = time.time()
    print "Loading genotypes from file took", end - start, "seconds."
    print "#Individuals:", len(genotypes), "   #SNPs:", len(genotypes[0])
    print "Memory cost of genotypes:", sys.getsizeof(genotypes)
    print
    
    print "---------- Phasing substring of first genotype ----------"
    lengthOfWindow = 120
    print "Window Size:", lengthOfWindow
    numOnes = 0
    for i in range(0, lengthOfWindow + 1):
        if genotypes[0][i] == 1:
            numOnes += 1
    print "Number of heterozygous positions:", numOnes
    start = time.time()
    phases = phaseGenotype(genotypes[0][:lengthOfWindow + 1])
    end = time.time()
    print "Number of phases:", len(phases)
    print "Memory cost of phases:", sys.getsizeof(phases)
    print "Took", end - start, "seconds."
    


        
        
