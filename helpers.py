from constants import aminoAcidMasses, params

def validTryptic(peptideSequence):
    return peptideSequence.endswith('K') or peptideSequence.endswith('R')

def endPeptide(peptideSequence):
    return validTryptic(peptideSequence)

def massOfPep(peptideSequence):
    return sum([aminoAcidMasses[c] for c in peptideSequence])

def goodPeptide(pep):
    mass, length = pep.neutralMass, len(pep.sequence)
    return mass >= params['MIN_PEPTIDE_MASS'] and mass <= params['MAX_PEPTIDE_MASS'] and length >= params['MIN_LEN_PEPTIDE'] and length <= params['MAX_LEN_PEPTIDE']


