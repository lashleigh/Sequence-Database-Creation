from constants import aminoAcidMasses, params

def checkForSpecialChar(char, numPhospho, numM):
    if char == 'P' or char == 'S' or char == 'T':
        numPhospho += 1
    elif char == 'M':
        numM += 1
    return(numPhospho, numM)

def badChar(c):
    return c == '*'

def validTryptic(peptideSequence):
    return peptideSequence.endswith('K') or peptideSequence.endswith('R')

def endPeptide(peptideSequence):
    return validTryptic(peptideSequence)

def massOfPep(peptideSequence):
    return sum([aminoAcidMasses[c] for c in peptideSequence])

def goodPeptide(neutralMass, sequence):
    mass, length = neutralMass, len(sequence)
    return mass >= params['MIN_PEPTIDE_MASS'] and mass <= params['MAX_PEPTIDE_MASS'] and length >= params['MIN_LEN_PEPTIDE'] and length <= params['MAX_LEN_PEPTIDE']

def enzymeChar(c):
    return c == 'K' or c == 'R'


