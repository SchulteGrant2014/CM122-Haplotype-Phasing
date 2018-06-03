#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 31 14:59:00 2018

@author: grant
"""

import time



# This function is the recursive helper function for phaseGenotype().
# It recursively computes all compatible haplotypes to the given genotype
# for all positions up to the current SNP, then phases the current SNP and
# updates the list of possible haplotypes to reflect including this new SNP.
# @param genotype The genotype to be phased
# @param i The current position in the genotype
# @return All compatible haplotype phases of genotype up to position i
def phaseGenotypeRecurse(genotype, i):
    
    # Base case: when we reach the beginning of the sequence
    if i == 0:
        if genotype[i] == 2:
            return [[1]]
        elif genotype[i] == 0:
            return [[0]]
        else:  # genotype[i] == 1
            return [[1], [0]]
    
    # Otherwise, get the phases leading up to our current SNP
    precedingPhases = phaseGenotypeRecurse(genotype, i-1)  # List of lists/haplotype phases
    
    # Depending on the copy number of the current SNP, either extend the phases (for 2/0) or add more haplotypes (for 1)
    currentSNP = genotype[i]
    phases = []
    if (currentSNP == 2):
        # Extend the current haplotype phases to include a 1 each
        for haplotype in range(0, len(precedingPhases)):
            newPhase = precedingPhases[haplotype]
            newPhase.append(1)
            phases.append(newPhase)
    
    elif (currentSNP == 0):
        # Extend the current haplotype phases to include a 0 each
        for haplotype in range(0, len(precedingPhases)):
            newPhase = precedingPhases[haplotype]
            newPhase.append(0)
            phases.append(newPhase)
    
    elif (currentSNP == 1):
        # Create 2 new haplotypes for each current haplotype
        for haplotype in range(0, len(precedingPhases)):
            phase1 = list(precedingPhases[haplotype])
            phase1.append(1)
            phase2 = precedingPhases[haplotype]
            phase2.append(0)
            phases.append(phase1)
            phases.append(phase2)
    else:
        raise ValueError("Genotype at index " + str(i) + " is not a valid symbol from the set {0, 1, 2}! Check that the provided genotype is correct. Error in phaseGenotypeRecurse().")
    
    # Return the newly constructed list of phases which includes the current SNP position, phased.
    return phases
    



# This function takes a genotype and outputs all compatible haplotype phases
# of that genotype.
# @param genotype The genotype to be phased
# @return All compatible haplotype phases of genotype
def phaseGenotype(genotype):
    # Get phases up to the last SNP
    lastIndex = len(genotype) - 1
    phases = phaseGenotypeRecurse(genotype, lastIndex)
    return phases
            



# This function takes a haplotype and a genotype and outputs the complementary
# haplotype to the one provided.
# @param haplotype One of the haplotypes making up the genotype
# @param genotype The genotype made up of the given haplotype and some
# other unknown haplotype
def haplotypeComplement(haplotype, genotype):
    
    # Make sure that the correct object type was provided to the function
    if not isinstance(haplotype, list):
        raise ValueError("Haplotype provided to haplotypeComplement() is not a list!")
    
    # Check that the haplotype and genotype provided are of the same length
    if len(haplotype) != len(genotype):
        raise ValueError("The lengths of Haplotype and Genotype arguments provided to haplotypeComplement() are not equal!")
    
    # Compare each haplotype position to the corresponding position in the genotype to derive the complementary haplotype
    complement = []
    for i in range(0, len(haplotype)):
        # Check that each position is compatible... if not, raise an error
        genotypeVal = genotype[i]
        haplotypeVal = haplotype[i]
        if ((genotypeVal == 2 and haplotypeVal != 1) or (genotypeVal == 0 and haplotypeVal != 0)):
            if not isinstance(haplotype[i], int):
                raise ValueError("Haplotype at SNP index " + str(i) + " is not of type int! Check that the provided haplotype is a single haplotype and not a list of all possible genotype phases. Error in haplotypeComplement().")
            else:
                raise ValueError("Haplotype is not compatible with genotype at index " + str(i) + " in function haplotypeComplement()!")
        else:
            complement.append(genotypeVal - haplotypeVal)
    
    # Return the complementary haplotype
    return complement




# Main function to run upon running this script
if __name__ == "__main__":
    
    
    # Run the genotype phasing algorithm
    start = time.time()
    
    #genotype = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    #genotype = [2, 1, 0, 1, 0, 1, 2, 2, 0]
    genotype = [2, 1, 1, 1, 0]
    phases = phaseGenotype(genotype)
    
    end = time.time()
    
    numOnes = 0
    for i in genotype:
        if i == 1:
            numOnes += 1
    
    print phases
    print "Total number of haplotypes:", len(phases)
    print "Length of genotype:", len(genotype)
    print "Number of heterozygous positions:", numOnes
    print "Phasing this genotype took " + str(end-start) + " seconds."
    
    """
    # Run the haplotype complement algorithm
    genotype = [2, 1, 0, 1, 0, 1, 2, 2, 0]
    haplotype = [1, 0, 0, 1, 0, 0, 1, 1, 0]
    complementaryHaplotype = haplotypeComplement(haplotype, genotype)
    
    print genotype
    print haplotype
    print complementaryHaplotype
    """




