#!/usr/bin/env python
from constants import aminoAcidMasses, params
from helpers import *
import sys
sys.setrecursionlimit(5500)

globalProteinList = []
globalPeptideList = []
globalGoodPeptides = []

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

def generateSemiTryptic(semiPeps):
    for pep in globalGoodPeptides:
        numPassedCleavages = 0
        numPassedPTS, numPassedM = 0, 0
        for i in range(1, len(pep.sequence)):
            if numPassedCleavages < max(1, pep.numRK - 1):
                if goodPeptide(pep.neutralMass, pep.sequence[i:]):
                    newPep = Peptide(pep.sequence[i:], massOfPep(pep.sequence[i:]), pep.numRK - numPassedCleavages, pep.numPTS - numPassedPTS, pep.numM - numPassedM, pep.peptideProteinList)
                    semiPeps.append(newPep)
                    #print numPassedCleavages, newPep.numRK, pep.sequence[:i], newPep
                    #print newPep.numPTS, newPep.numM, pep.sequence[:i], pep.sequence[i:]
                if enzymeChar(pep.sequence[i]):
                    numPassedCleavages += 1
                numPassedPTS, numPassedM = checkForSpecialChar(pep.sequence[i], numPassedPTS, numPassedM)

def findNextPeptide(proteinSequence, protein):
    pepSeq, numPTS, numM = '', 0, 0
    for i, char in enumerate(proteinSequence):
        if badChar(char):
            continue
        pepSeq += char
        numPTS, numM = checkForSpecialChar(char, numPTS, numM)
        if endPeptide(pepSeq):
            newPep = Peptide(pepSeq, massOfPep(pepSeq), 1, numPTS, numM, [protein])
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
                if goodPeptide(potentialPep.neutralMass, potentialPep.sequence):
                    proteinPeptideList.append(potentialPep)
        proteinPeptideList += findGoodPeptides(temp[1:])
    return proteinPeptideList 

class Protein(object):
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence
        globalProteinList.append(self)

    def __repr__(self):
        return self.name[:7]

class Peptide(object):
    def __init__(self, sequence = '', neutralMass = 0, numRK = 0, numPTS = 0, numM = 0, peptideProteinList = set()):
        self.sequence = sequence
        self.neutralMass = neutralMass
        self.numRK = numRK
        self.numPTS = numPTS
        self.numM = numM 
        self.peptideProteinList = peptideProteinList

    def __add__(self, other):
        sequence = self.sequence + other.sequence
        neutralMass = self.neutralMass + other.neutralMass
        numRK = self.numRK + other.numRK
        numPTS = self.numPTS + other.numPTS
        numM = self.numM + other.numM 
        proteins = self.peptideProteinList.union(other.peptideProteinList)
        return Peptide(sequence, neutralMass, numRK, numPTS, numM, proteins)

    def __repr__(self):
        return self.sequence

    def __eq__(self, other):
        return self.sequence == other.sequence and \
                self.neutralMass == other.neutralMass and \
                self.numRK == other.numRK and \
                self.numPTS == other.numPTS and \
                self.numM == other.numM

    def merge(self, other):
        self.peptideProteinList.extend(other.peptideProteinList)


if len(sys.argv) < 2:
    print "\n\t###  Please supply one or more .fasta files for digestion"
    print "\t###  example:   ./peptideCreator.py test.fasta\n"
    sys.exit(1)

for fname in sys.argv[1:]:
    for protein in get_proteins(fname):
        temp = digest(protein.sequence, protein)
        #print protein.sequence
        protein.peptides = findGoodPeptides(temp)
        globalGoodPeptides.extend(protein.peptides)

if params['SEMI_TRYPTIC']:
    semiPeps = []
    generateSemiTryptic(semiPeps)
    globalGoodPeptides.extend(semiPeps)
#for pep in sorted(globalGoodPeptides, key = lambda peptide: peptide.neutralMass):
#    print pep.numPTS, pep.numM, pep.numRK, pep
    
#for pep in globalGoodPeptides:
#    print pep.numPTS, pep.numM, pep.numRK, pep
