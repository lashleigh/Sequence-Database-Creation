#!/usr/bin/env python
from constants import aminoAcidMasses, params
from helpers import *
import sys

globalProteinList = []
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

def badChar(c):
    return c == '*'

def findNextPeptide(proteinSequence, protein):
    pepSeq, numPTS, numM = '', 0, 0
    for i, char in enumerate(proteinSequence):
        if badChar(char):
            continue
        pepSeq += char
        numPTS += checkPhoso(char)
        numM += checkMeth(char)
        if endPeptide(pepSeq):
            newPep = Peptide(pepSeq, massOfPep(pepSeq), numPTS, numM, [protein])
            globalPeptideList.append(newPep)
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
                    proteinPeptideList.append(potentialPep)
        proteinPeptideList += findGoodPeptides(temp[1:])
    return proteinPeptideList 

class Protein(object):
    def __init__(self, name, sequence, peptides = []):
        self.name = name
        self.sequence = sequence
        self.peptides = []
        globalProteinList.append(self)

    def __repr__(self):
        return self.name[:7]

class Peptide(object):
    def __init__(self, sequence = '', neutralMass = 0, numPTS = 0, numM = 0, peptideProteinList = set()):
        self.sequence = sequence
        self.neutralMass = neutralMass
        self.numPTS = numPTS
        self.numM = numM 
        self.peptideProteinList = peptideProteinList

    def __add__(self, other):
        sequence = self.sequence + other.sequence
        neutralMass = self.neutralMass + other.neutralMass
        numPTS = self.numPTS + other.numPTS
        numM = self.numM + other.numM 
        proteins = self.peptideProteinList.union(other.peptideProteinList)
        return Peptide(sequence, neutralMass, numPTS, numM, proteins)

    def __repr__(self):
        return self.sequence

if len(sys.argv) < 2:
    print "\n\t###  Please supply one or more .fasta files for digestion"
    print "\t###  example:   ./peptideCreator.py test.fasta\n"
    sys.exit(1)

for fname in sys.argv[1:]:
    for protein in get_proteins(fname):
        temp = digest(protein.sequence, protein)
        #print protein.sequence
        protein.peptides = findGoodPeptides(temp)
        for pep in sorted(protein.peptides, key = lambda peptide: peptide.neutralMass):
            print pep.neutralMass,'\t', pep.numPTS, pep.numM, pep, pep.peptideProteinList
    




