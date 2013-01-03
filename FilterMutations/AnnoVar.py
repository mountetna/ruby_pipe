#!/usr/bin/env python

import os
import subprocess
import glob

AnnoVarSuffix = ".variant_function"
AnnoVarExonicSuffix = ".exonic_variant_function"
VarLineLen = 7
ExonicLineLen = 8
NA = "NA"

class AnnoVarLine:
	''' This class is used to hold an AnnoVar element and print it out properly '''
	def __init__(self, varLine, annoVarFilters, exonicLine=None):
		''' Initialize an AnnoVarLine instance, incorporating exonic mutation information if applicable. '''
		if varLine.endswith('\n'):
			varLine = varLine[:-1]
		varLine = varLine.split('\t')
		if len(varLine) != VarLineLen:
			errAbort("Expecting AnnoVar %s line to have %d fields, not %d: %s" % (AnnoVarSuffix, VarLineLen, len(varLine), "  ".join(varLine)))
		self.context = varLine[0]
		self.genes = self.parseGenes(varLine[1])
		self.chrom = varLine[2]
		self.startPos = int(varLine[3])
		self.endPos = int(varLine[4])
		self.ref_allele = varLine[5]
		self.alt_allele = varLine[6]
		self.knownVariant = "NOVEL"
		self.nucleotide = NA
		self.protein = NA
		self.accession = NA
		self.exon = NA
		self.exonMutType = NA

		if exonicLine != None:
			self.nucleotide = exonicLine.nucleotides
			self.protein = exonicLine.aminoacids
			self.accession = exonicLine.accessions
			self.exon = exonicLine.exons
			self.exonMutType = exonicLine.exonMutType

		self.makeOurJudgment(annoVarFilters)

	def parseGenes(self, s):
		''' Parse the "genes" column of the .variant_function file to be something we want to use '''
		# Uncomment this to have the previously-used output
		# return s
		
		# Find all the open- and close-parentheses, so we can omit all ',' characters within them
		opIdx = [] ; cpIdx = []
		for idx, c in enumerate(s):
			if c == '(':
				opIdx.append(idx)
			elif c == ')':
				cpIdx.append(idx)

		numParens = len(opIdx)
		if numParens == 0:
			return s.replace(";", ",")

		# Verify parens are well-formatted and not nested.
		assert numParens == len(cpIdx)
		for i in range(numParens):
			assert opIdx[i] < cpIdx[i]
			if i < numParens-1:
				assert cpIdx[i] < opIdx[i+1]

		start = 0
		retval = ""
		for i in range(numParens):
			retval += s[start:opIdx[i]]
			start = cpIdx[i]+1
		retval = ",".join([x for x in retval.split(',') if x != "NONE"])
		return retval.replace(";", ",")

	def getKnownVariantStatus(self):
		return self.knownVariant

	def makeOurJudgment(self, annoVarFilters):
		self.ourJudgment = "yes"
		judgmentReasons = []
		numContexts = 0
		numBadContexts = 0
		for context in self.context.split(";"):
			numContexts += 1
			if context in annoVarFilters['contextsToOmit']:
				numBadContexts += 1

		if numBadContexts == numContexts:
			self.ourJudgment = "no"
			judgmentReasons.append("ANNOTATION_CONTEXT_TO_OMIT")

		self.ourJudgmentReasons = ",".join(judgmentReasons)

	def shouldBeKept(self, annoVarRemovalFilters):
		''' Return True iff we should report this alteration '''
		return True

	def exonMutTypeStr(self):
		''' Return the exonic mutation type '''
		if self.exonMutType == NA:
			return NA
		else:
			if self.exonMutType.startswith('synonymous'):
				return 'Silent'
			elif self.exonMutType.startswith('stopgain'):
				return 'Nonsense'
			elif self.exonMutType.startswith('stoploss'):
				return 'Non-stop'
			elif self.exonMutType.startswith('nonsynonymous'):
				return 'Missense'
			elif self.context.startswith('splicing'):
				return 'Splice-site'
			return self.exonMutType

class ExonicAnnoVarLine:
	def __init__(self, line):
		if line.endswith('\n'):
			line = line[:-1]
		line = line.split('\t')
		if len(line) != ExonicLineLen:
			errAbort("Expecting AnnoVar %s line to have %d fields, not %d: %s" % (AnnoVarExonicSuffix, ExonicLineLen, len(line), "  ".join(line)))
		assert line[0].startswith('line')
		self.lineNum = int(line[0].replace('line', '', 1))
		self.exonMutType = line[1]
		self.chrom = line[3]
		self.startPos = line[4]
		self.endPos = line[5]
		self.ref_allele = line[6]
		self.alt_allele = line[7]
	
		# Parse the 3rd column for all of its information.  This is a set of comma-delimited entries (one for each gene affected),
		# that has the following 5 columns delimited by ':':  gene, accession, exon #, nucleotide changed (with position in gene),
		# amino acid changed (with position in gene)
		unk = 'UNKNOWN'
		if line[2] == 'UNKNOWN':
			self.genes = unk
			self.accessions = unk
			self.exons = unk
			self.nucleotides = unk
			self.aminoacids = unk
		else:
			entries = [x.split(':') for x in line[2].split(',') if len(x) > 0]
			self.genes = ','.join([x[0] for x in entries])
			self.accessions = ','.join([x[1] for x in entries])
			def getExon(s):
				if s == 'wholegene':
					return s
				return int(s.replace('exon','',1))
			self.exons = ','.join([str(getExon(x[2])) for x in entries])
			self.nucleotides = ','.join([x[3].replace('c.','',1) for x in entries if len(x) > 3])
			if self.nucleotides == "":
				self.nucleotides = NA
			self.aminoacids = ','.join([x[4].replace('p.','',1) for x in entries if len(x) > 4])
			if self.aminoacids == "":
				self.aminoacids = NA

def getAnnoVarResults(annoVarInputFn, genome, executable, annoVarDBDir, thousandGenomesDBs, dbSNPDBs, annoVarFilters):
	''' Return AnnoVar objects for the results '''
	# First run the raw AnnoVar script to get all variations
	cmd = executable + " --buildver %s --geneanno %s %s" % (genome, annoVarInputFn, annoVarDBDir)
	devnull = open(os.devnull)
	proc = subprocess.Popen(cmd, shell=True, stderr=devnull)
	proc.communicate()

	retval = []
	allVarFile = open(annoVarInputFn + AnnoVarSuffix)
	exonVarFile = open(annoVarInputFn + AnnoVarExonicSuffix)
	nextExonVarStr = exonVarFile.readline()
	if not nextExonVarStr:
		nextExonVar = None
		exonVarFile.close()
	else:
		nextExonVar = ExonicAnnoVarLine(nextExonVarStr)
	allVarLineNum = 0
	prevExonLineNum = -1

	print "Screening AnnoVar outputs with the following parameters:"
	print "  Annotating AnnoVar outputs with the following:"
	for k in sorted(annoVarFilters.keys()):
		print "    ", k, ":", annoVarFilters[k]

	annoVarFilters['contextsToOmit'] = set(annoVarFilters['contextsToOmit'].split(','))

	# Loop through the entire .variant_function file, making an entry for each.  For those that are exonic, load
	# the proper exonic variant data as well
	while True:
		varLine = allVarFile.readline()
		allVarLineNum += 1
		if not varLine:
			break
		if (nextExonVar != None) and (nextExonVar.lineNum == allVarLineNum):
			a = AnnoVarLine(varLine, annoVarFilters, nextExonVar)

			assert nextExonVar.lineNum > prevExonLineNum
			prevExonLineNum = nextExonVar.lineNum

			nextExonVarStr = exonVarFile.readline()
			if not nextExonVarStr:
				nextExonVar = None
				exonVarFile.close()
			else:
				nextExonVar = ExonicAnnoVarLine(nextExonVarStr)
		else:
			a = AnnoVarLine(varLine, annoVarFilters, None)

		retval.append(a)
			
	if nextExonVar != None:
		errAbort("All exonic variants should have been associated with their variant already.")
	allVarFile.close()
	

	# Now run the 1000 genomes- and dbSNP-filtered runs to see which mutations are known germline mutations
	for db in [x for x in dbSNPDBs.split(',') + thousandGenomesDBs.split(',') if len(x) > 0]:
		cmd = executable + " --buildver %s -filter -dbtype %s %s %s" % (genome, db, annoVarInputFn, annoVarDBDir)
		print "working on", db
		proc = subprocess.Popen(cmd, shell=True, stderr=devnull)
		proc.communicate()
	devnull.close()

	# Build up a set of mutations that are found in either dbSNP or 1000 genomes databases that can be omitted
	knownVariants = {}
	for inNormalVariantsFn in glob.glob(annoVarInputFn + "*_dropped"):
		for line in open(inNormalVariantsFn):
			db, name, chrom, start, end, ref, alt = line.strip().split('\t')
			knownVariants["%s:%s-%s" % (chrom, start, end)] = knownVariants.get("%s:%s-%s" % (chrom, start, end), []) + [db + "_" + name]

	for variant in retval:
		v = knownVariants.get("%s:%d-%d" % (variant.chrom, variant.startPos, variant.endPos), [])
		if v != []:
			v.sort()
			variant.knownVariant = ",".join(v)

	return retval
