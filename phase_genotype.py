#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 31 14:59:00 2018

@author: grant
"""


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
    
    else:  # if currentSNP == 1
        # Create 2 new haplotypes for each current haplotype
        for haplotype in range(0, len(precedingPhases)):
            phase1 = list(precedingPhases[haplotype])
            phase1.append(1)
            phase2 = precedingPhases[haplotype]
            phase2.append(0)
            phases.append(phase1)
            phases.append(phase2)
    
    # Return the newly constructed list of phases which includes the current SNP position, phased.
    return phases
    


def phaseGenotype(genotype):
    # Get phases up to the last SNP
    lastIndex = len(genotype) - 1
    phases = phaseGenotypeRecurse(genotype, lastIndex)
    return phases
            


"""
if __name__ == "__main__":
    genotype = [2, 0, 1, 1, 2, 2, 1]
    print phaseGenotype(genotype)
"""



