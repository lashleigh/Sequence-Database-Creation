#!/usr/bin/env python
from constants import aminoAcidMasses, params
from helpers import *
import sys

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

def findNextPeptide(seq):
    pep = ''
    for i, char in enumerate(seq):
        pep += char
        if endPeptide(pep):
            return(pep, seq[i+1:])
    return(None, None)

def digest(sequence):
    peptides = []
    if sequence:
        peptide, next_sequence = findNextPeptide(sequence)
        if peptide:
            peptides.append(peptide)
            if next_sequence:
                peptides += digest(next_sequence)
    return peptides 

def findGoodPeptides(peptides, temp):
    if temp:
        potentialPep = ''
        for i in range(params['ALLOWED_MISSED_CLEAVAGES'] + 1):
            if i < len(temp):
                potentialPep += temp[i]
                if goodPeptide(potentialPep):
                    peptides.append(potentialPep)
                    print i, potentialPep
            else:
                break
        findGoodPeptides(peptides, temp[1:])

class Protein(object):
    def __init__(self, name, sequence, peptides = []):
        self.name = name
        self.sequence = sequence
        self.peptides = []
        proteins.append(self)

class Peptide(object):
    def __intit__(self, protein, sequence, neutralMass = 0, numPTS = 0, numM = 0):
        self.protein = []
        self.sequence = sequence

if len(sys.argv) < 2:
    print "\n\t###  Please supply one or more .fasta files for digestion"
    print "\t###  example:   ./peptideCreator.py test.fasta\n"
    sys.exit(1)

for fname in sys.argv[1:]:
    for p in get_proteins(fname):
        temp = digest(p.sequence)
        findGoodPeptides(p.peptides, temp)
        #if params['ALLOWED_MISSED_CLEAVAGES']:
        #    p.peptides = miscleave(p.peptides, params['ALLOWED_MISSED_CLEAVAGES'])
        print sorted(p.peptides)
        print
        #for pep in p.peptides:
        #    print pep, massOfPep(pep)

