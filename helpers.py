from constants import aminoAcidMasses, params

def validTryptic(pep):
    return pep.endswith('K') or pep.endswith('R')

def endPeptide(pep):
    return validTryptic(pep)

def massOfPep(peptide):
    return sum([aminoAcidMasses[c] for c in peptide])

def goodPeptide(pep):
    mass, length = massOfPep(pep), len(pep)
    return mass >= params['MIN_PEPTIDE_MASS'] and mass <= params['MAX_PEPTIDE_MASS'] and length >= params['MIN_LEN_PEPTIDE'] and length <= params['MAX_LEN_PEPTIDE']


