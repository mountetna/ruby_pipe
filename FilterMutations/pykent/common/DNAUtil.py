#!/usr/bin/env python

"""
A utility class for functions on DNA
"""

from Sanity import errAbort

DNA = set('ACGTacgtNn-')

dnaRC = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'a':'t', 't':'a', 'g':'c', 'c':'g', 'N':'N', 'n':'n'}

codonTable = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', \
			  'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', \
			  'TAT':'Y', 'TAC':'Y', 'TAA':0, 'TAG':0, \
			  'TGT':'C', 'TGC':'C', 'TGA':0, 'TGG':'W', \
\
			  'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', \
			  'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', \
			  'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q', \
			  'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', \
\
			  'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', \
			  'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T', \
			  'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K', \
			  'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R', \
\
			  'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V', \
			  'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A', \
			  'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', \
			  'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G', \
}

NON_CODON = 'X'
SILENT_MUT = "Silent"
MISSENSE_MUT = "Missense"
NONSENSE_MUT = "Nonsense"
GAIN_MUT = "Gain Mutation"
FS_INS = "Frameshift Insertion"
INFRAME_INS = "Inframe Insertion"
FS_DEL = "Frameshift Deletion"
INFRAME_DEL = "Inframe Deletion"

def isDNA(s):
	''' Return True iff the input string is entirely composed of DNA-like characters including N or '-' '''
	for c in s:
		if c not in DNA:
			return False
	return True

def complement(dna):
	''' Return the complement DNA (not reversed) '''
	return "".join([dnaRC.get(x,x) for x in dna])

def reverseComplement(dna):
	''' Return the reverse complement DNA '''
	return "".join([dnaRC.get(x,x) for x in reversed(dna)])

def toRNA(dna):
	''' Return an RNA sequence from DNA (convert T to U) '''
	return dna.replace('t', 'u').replace('T','U')

def codonValue(dna):
	'''Return the codon value of the dna string or 'X' if it is not a valid DNA string '''
	global codonTable, NON_CODON
	if len(dna) != 3:
		errAbort("codonValue takes a 3-bp dna string as argument: %s" % dna)

	return codonTable.get(dna.upper(), NON_CODON)

def isStopCodon(dna):
	''' Return True iff dna is a stop codon '''
	return codonValue(dna) == 0

def exonMutationType(wtDNA, mutantDNA):
	''' Return the type of mutation.  Assumes both wt and mutant are fully within exons '''
	global NON_CODON, SILENT_MUT, MISSENSE_MUT, NONSENSE_MUT, GAIN_MUT, FS_INS, INFRAME_INS, FS_DEL, INFRAME_DEL
	wtLen = len(wtDNA.replace('-',''))
	mutantLen = len(mutantDNA.replace('-',''))
	if wtLen < mutantLen:
		if (mutantLen-wtLen)%3 != 0:
			return FS_INS
		else:
			return INFRAME_INS
	elif wtLen > mutantLen:
		if (wtLen-mutantLen)%3 != 0:
			return FS_DEL
		else:
			return INFRAME_DEL
	else:
		wtCodon = codonValue(wtDNA)
		mutantCodon = codonValue(mutantDNA)
		if (wtCodon == NON_CODON) or (mutantCodon == NON_CODON):
			errAbort("Invalid DNA codons tested for mutation type.  Wt: %s, mutant: %s" % (wtDNA, mutantDNA))

		if not isStopCodon(wtDNA):
			if isStopCodon(mutantDNA):
				return NONSENSE_MUT
			else:
				if wtCodon == mutantCodon:
					return SILENT_MUT
				else:
					return MISSENSE_MUT
		else:
			if isStopCodon(mutantDNA):
				return SILENT_MUT
			else:
				return GAIN_MUT
