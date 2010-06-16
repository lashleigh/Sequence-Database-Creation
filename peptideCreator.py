#!/usr/bin/env python
import math
import sys
from constants import aminoAcidMasses, params
from pprint import pprint

proteins = []

def get_proteins(fname):
    name, sequence = None, None
    for line in open(fname):
        if line.startswith(">"):
            if name:
                yield Protein(name, sequence)
            name = line
            sequence = ""
        else:
            sequence += line.strip()
    yield Protein(name, sequence)

def validTryptic(pep):
    return pep.endswith('K') or pep.endswith('R')

def endPeptide(pep):
    return validTryptic(pep)

def massOfPep(peptide):
    return sum([aminoAcidMasses[c] for c in peptide])

def goodPeptide(pep):
    mass, length = massOfPep(pep), len(pep)
    return mass >= params['MIN_PEPTIDE_MASS'] and mass <= params['MAX_PEPTIDE_MASS'] and length >= params['MIN_LEN_PEPTIDE'] and length <= params['MAX_LEN_PEPTIDE']

def findNextPeptide(seq):
    pep = ''
    for i, char in enumerate(seq):
        pep += char
        if endPeptide(pep):
            return(pep, seq[i+1:])
    return(None, None)

def digest(sequence, peptides):
    if sequence:
        peptide, next_sequence = findNextPeptide(sequence)
        if peptide:
            if goodPeptide(peptide):
                peptides.append(peptide)
            if next_sequence:
                digest(next_sequence, peptides)

class Protein(object):
    def __init__(self, name, sequence, peptides = []):
        self.name = name
        self.sequence = sequence
        self.peptides = []
        proteins.append(self)

for fname in sys.argv[1:]:
    for p in get_proteins(fname):
        digest(p.sequence, p.peptides)
        #if params['ALLOWED_MISSED_CLEAVAGES']:
        #    p.peptides = miscleave(p.peptides, params['ALLOWED_MISSED_CLEAVAGES'])
        print p.peptides
        #for pep in p.peptides:
        #    print pep, massOfPep(pep)

