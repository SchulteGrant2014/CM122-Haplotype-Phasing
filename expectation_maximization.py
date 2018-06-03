#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  2 14:52:04 2018

@author: grant
"""

from phase_genotype import phaseGenotype
from phase_genotype import haplotypeComplement
from load_genome import load_genome_faster
import time


def calcHaplotypeProbability(expectation, numGenotypes):
    return expectation / (2 * numGenotypes)



def calcHaplotypePairingExpectation(haplotypePairProbability, totalProbability):
    return haplotypePairProbability / totalProbability



def calcHaplotypePairingProbability(haplotype_1, haplotype_2, haplotypeProbabilityDict):
    haplotype_1_str = ''.join(map(str, haplotype_1))
    haplotype_2_str = ''.join(map(str, haplotype_2))
    haplo_1_prob = haplotypeProbabilityDict[haplotype_1_str]
    haplo_2_prob = haplotypeProbabilityDict[haplotype_2_str]
    return haplo_1_prob * haplo_2_prob  # Probability of two haplotypes forming the genotype is the frequency of both of the haplotypes in the population multiplied (probability of having 1 AND 2)




def expectation_maximization(genotypes, numberOfIterations):
    
    # Error checking... make sure "genotypes" is of correct type: list of multiple genotypes = list of lists of SNPs
    if isinstance(genotypes, list):
        if len(genotypes) == 0:
            return []
        if isinstance(genotypes[0], list):
            if not isinstance(genotypes[0][0], int):
                raise ValueError("Genotype list provided to expectation_maximization() is not of the correct type! Type should be list of genotypes, which is a list of lists (of SNP=int).")
        else:
            raise ValueError("Genotype list provided to expectation_maximization() is not of the correct type! Type should be a list of many genotypes, which is a list of lists (of SNP=int).")
    else:
        raise ValueError("Genotype list provided to expectation_maximization() is not of the correct type! Type should be list of genotypes, which is a list of lists (of SNP=int).")
    
    
    # Phase the genotypes into their possible phases/haplotypes.
    # Use complements to see if any pairings are "duplicates"
    # i.e. order of the haplotypes doesn't matter: [0,1] and [1,0] is equivalent to [1,0] and [0,1]
    genotypePhases = []  # For each genotype, a list of every possible phase WITHOUT "duplicates"
    for geno in genotypes:
        possibleTopHaplotypes = []  # The set of all possible phases, excluding those whose self/complement ordering is merely the reverse of another one already in the set
        discoveredHaploDict = dict()  # Keep track of which haplotype pairs have already been added to the list of phases for this genotype... ordering of haplotypes is irrelevant: don't add [[0,1],[1,0]] if [[1,0],[0,1]] is already accounted for!
        phases = phaseGenotype(geno)
        for phase in phases:  # Add each haplotype to the list of possible phases for this genotype ONLY IF it (or its complement) isn't already discovered
            complement = haplotypeComplement(phase, geno)
            # Make sure we don't duplicate entries - order of haplotypesdoesn't matter
            phase_str = ''.join(map(str, phase))
            complement_str = ''.join(map(str, complement))
            if phase_str not in discoveredHaploDict:
                discoveredHaploDict[phase_str] = 1
                discoveredHaploDict[complement_str] = 1
                possibleTopHaplotypes.append(phase)  # <===============
                
        genotypePhases.append(possibleTopHaplotypes)  # WITHOUT "duplicates"
    # Note: genotypePhases[i] = the i_th genotype; genotypePhases[i][j] = one compatible top haplotype of the i_th genotype, which is never the same as another haplotype's complement
    
    
    # -------------------- Set up the haplotype probability table --------------------
    
    # First, count how many unique haplotypes exist from the phasing
    haplotypeProbDict = dict()  # Initialize a dictionary of all haplotypes and their probabilities
    numHaplotypes = 0
    for geno in range(0, len(genotypePhases)):
        for haplo in genotypePhases[geno]:
            haplo_str = ''.join(map(str, haplo))  # Convert the haplotype (list) to a string
            if haplo_str not in haplotypeProbDict:
                numHaplotypes += 1
                haplotypeProbDict[haplo_str] = 0.0  # Add haplotype to dict to mark it as "discovered"
            complement_str = ''.join(map(str, haplotypeComplement(haplo, genotypes[geno])))
            if complement_str not in haplotypeProbDict:
                numHaplotypes += 1
                haplotypeProbDict[complement_str] = 0.0  # Add haplotype to dict to mark it as "discovered"
    
    # Update the probabilities of each haplotype to be equal, initially (our "guess")
    initial_haplo_prob = (1.0 / numHaplotypes)
    for haplo in haplotypeProbDict.iterkeys():
        haplotypeProbDict[haplo] = initial_haplo_prob  # Initialize haplotype probability to [1 / (number of unique haplotypes)]
    
    
    # -------------------- EM (Expectation/Maximization Algorithm)  --------------------
    
    # Set up the haplotype pair probabilities (middle column from slides) to be the correct dimensions and have initial parameters
    haplotypePairProbs = []
    for geno in genotypePhases:
        numPhases = len(geno)
        pairProbabilities = [ (1.0 / numPhases) for t in range(0, numPhases) ]
        haplotypePairProbs.append(pairProbabilities)
    
    
    # Update the right column (haplo probs), then the middle (pair probs) = 1 iteration
    for iteration in range(1, numberOfIterations+1):
        
        print "--------------- Iteration:", iteration, "---------------"
        # Update the right column (haplotype frequencies)
        for haplotypeOfInterest in haplotypeProbDict.iterkeys():
            haplo_list = list(map(int, haplotypeOfInterest))  # Convert back to list (haplotype format) from string (key format)
            sumOfPairProbs = 0.0
            for geno in range(0, len(genotypePhases)):  # Look through each individual to see everywhere the haplotype shows up, keep rolling sum of probabilities
                for haplo in range(0, len(genotypePhases[geno])):  # Check each haplotype from the individual to see if it or its complement is here
                    if (genotypePhases[geno][haplo] == haplo_list) or (haplotypeComplement(genotypePhases[geno][haplo], genotypes[geno]) == haplo_list):
                        sumOfPairProbs += haplotypePairProbs[geno][haplo]
            
            # Update the right column (haplotype pair probabilities) given the newly-calculated expectation
            haplotypeProbDict[haplotypeOfInterest] = calcHaplotypeProbability(sumOfPairProbs, len(genotypes))
            
        # Print the results to the console of updating the haplotype probabilities
        for haplotype in haplotypeProbDict.iterkeys():
            print haplotype + ":", haplotypeProbDict[haplotype]
        
        # Update the middle column using the newly-computed right-column
        for geno in range(0, len(genotypePhases)):
            pairingProbs = [ 0.0 for t in range(0, len(genotypePhases[geno]))]
            totalProb = 0.0
            for haplo in range(0, len(genotypePhases[geno])):
                haplotype = genotypePhases[geno][haplo]
                genotype = genotypes[geno]
                complement = haplotypeComplement(haplotype, genotype)
                pairingProbs[haplo] = calcHaplotypePairingProbability(haplotype, complement, haplotypeProbDict)
                totalProb += pairingProbs[haplo]
            
            # After calculating all individual pairing probabilities, calculate the expectation of each
            for haplo in range(0, len(genotypePhases[geno])):
                haplotypePairProbs[geno][haplo] = calcHaplotypePairingExpectation(pairingProbs[haplo], totalProb)
    
    
        print "Haplotype Pair Probabilities:", haplotypePairProbs
    
    
    # At the end, go through all of the possible phases and choose the maximum for each genotype
    for genotype in range(0, len(genotypePhases)):
        maxProbabilityPhase = []
        maxProbability = 0.0
        for phase in range(0, len(genotypePhases[genotype])):
            phase_list = genotypePhases[genotype][phase]
            phase_str = ''.join(map(str, phase_list))
            probability = haplotypePairProbs[genotype][phase]
            if probability > maxProbability:
                maxProbability = probability
                maxProbabilityPhase = phase_list
        maxProbabilityPhaseComplement = haplotypeComplement(maxProbabilityPhase, genotypes[genotype])
        print "Maximum probability phase for Genotype", genotype, "=", maxProbability
        print maxProbabilityPhase
        print maxProbabilityPhaseComplement



if __name__ == "__main__":
    
    start = time.time()
    genos = [[2, 1, 1, 1, 0], [1, 0, 0, 0, 1], [2, 2, 2, 2, 1]]  # Example from lecture "EM Final", slide 41
    numberOfIterations = 20
    expectation_maximization(genos, numberOfIterations)
    end = time.time()
    
    print "This algorithm took", end - start, "seconds to run."
    
    
    
    









