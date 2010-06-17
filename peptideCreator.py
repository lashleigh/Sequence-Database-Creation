#!/usr/bin/env python
from constants import aminoAcidMasses, params
from helpers import *
import sys

proteins = []
globalPeptideList = []

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

def findNextPeptide(proteinSequence, protein):
    pepSeq, numPTS, numM = '', 0, 0
    for i, char in enumerate(proteinSequence):
        pepSeq += char
        numPTS += checkPhoso(char)
        numM += checkMeth(char)
        if endPeptide(pepSeq):
            newPep = Peptide(protein, pepSeq, massOfPep(pepSeq), numPTS, numM)
            return(newPep, proteinSequence[i+1:])
    return(None, None)

def digest(proteinSequence, protein):
    tempPeptideList = []
    if proteinSequence:
        peptide, next_sequence = findNextPeptide(proteinSequence, protein)
        if peptide:
            tempPeptideList.append(peptide)
            if next_sequence:
                tempPeptideList += digest(next_sequence, protein)
    return tempPeptideList 

def findGoodPeptides(temp):
    proteinPeptideList = []
    if temp:
        potentialPep = Peptide()
        for i in range(params['ALLOWED_MISSED_CLEAVAGES'] + 1):
            if i < len(temp):
                potentialPep += temp[i]
                if goodPeptide(potentialPep):
                    print i, potentialPep.sequence
                    proteinPeptideList.append(potentialPep)
            else:
                break
        proteinPeptideList += findGoodPeptides(temp[1:])
    return proteinPeptideList 

class Protein(object):
    def __init__(self, name, sequence, peptides = []):
        self.name = name
        self.sequence = sequence
        self.peptides = []
        proteins.append(self)

class Peptide(object):
    def __init__(self, protein = '', sequence = '', neutralMass = 0, numPTS = 0, numM = 0):
        self.protein = protein
        self.sequence = sequence
        self.neutralMass = neutralMass
        self.numPTS = numPTS
        self.numM = numM 
        globalPeptideList.append(self)

    def __add__(self, other):
        self.protein = other.protein
        self.sequence += other.sequence
        self.neutralMass += other.neutralMass
        self.numPTS += other.numPTS
        self.numM += other.numM 
        return Peptide(self.protein, self.sequence, self.neutralMass, self.numPTS, self.numM)


if len(sys.argv) < 2:
    print "\n\t###  Please supply one or more .fasta files for digestion"
    print "\t###  example:   ./peptideCreator.py test.fasta\n"
    sys.exit(1)

for fname in sys.argv[1:]:
    for protein in get_proteins(fname):
        temp = digest(protein.sequence, protein)
        protein.peptides = findGoodPeptides(temp)
        for pep in protein.peptides:
            print pep.sequence, pep.neutralMass
