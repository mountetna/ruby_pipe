#!/usr/bin/env python

"""
Filter.py - Restrict input mutation calls to ones we want to keep.

Usage: Filter.py mutationConfig.cfg pointMutationFile.txt indelFile.txt output.txt
 where
 mutationConfig.cfg     holds the configuration for all types of mutation callers, tool locations, etc.
 pointMutationFile.txt  holds a list of all point mutations identified
 indelFile.txt          holds a list of all insertion/deletion alterations identified
 output.txt             holds the output data with all filters applied

Options:
  -h, -?, --help       Print this help and exit
  --keepTmpFiles       If flag set, do not delete intermediate AnnoVar output files (default=False)
  --tmp <dir>          Use the given directory, otherwise use the current working directory
"""

import sys
import os
import getopt
import ConfigParser
import tempfile
import datetime
import glob
from pykent.common.Sanity import errAbort, canBeInt, canBeNum, roughlyEqual
from pykent.common.DNAUtil import isDNA
from MuTector import *
from SomaticIndel import *
from AnnoVar import *

# Global variables
NUM_ARGS = 4
KEEP_TMP_FILES = False
TMP_DIR="."

BIG_INT = 99999999
BIG_FLOAT = 99999999.9

# All screening filters
AllFilters = {'SNPRemovalFilters': {'muTectJudgmentsAllowed':"KEEP"}, \
              'SNPAnnotationFilters': {'minTumorReads':0, 'minTumorAltReads':0, 'minTumorVarFreq':0.0, 'minNormalReads':0, 'maxNormalAltReads':BIG_INT, 'maxNormalVarFreq':1.0, 'coveredAllowed':"", 'normalVariantFreq': '>= tumorVariantFreq'}, \
              'IndelRemovalFilters': {'minTumorReads':0, 'minTumorAltReads':0, 'minTumorVarFreq':0.0, 'minNormalReads':0, 'maxNormalAltReads':BIG_INT, 'maxNormalVarFreq':1.0, \
                'maxAvgTumorAltMismatch':BIG_FLOAT, 'maxAvgTumorRefMismatch':BIG_FLOAT, 'maxAvgNormalAltMismatch':BIG_FLOAT, 'maxAvgNormalRefMismatch':BIG_FLOAT, 'minAvgTumorAltBaseQualInNQS':0.0, 'minAvgTumorRefBaseQualInNQS':0.0, \
                'minAvgNormalAltBaseQualInNQS':0.0, 'minAvgNormalRefBaseQualInNQS':0.0, 'maxAvgTumorAltMismatchInNQS':BIG_FLOAT, 'maxAvgTumorRefMismatchInNQS':BIG_FLOAT, 'maxAvgNormalAltMismatchInNQS':BIG_FLOAT, 'maxAvgNormalRefMismatchInNQS':BIG_FLOAT, \
                'normalVariantFreq': '>= tumorVariantFreq'}, \
			  'IndelAnnotationFilters': {'maxHomopolymerRuns':BIG_INT}, \
              'AnnoVarRemovalFilters': {}, \
	      'AnnoVarAnnotationFilters': {'contextsToOmit':"",'mutTypesToOmit':""}}

def configure(fn):
	''' Read the configuration file '''
	global AllFilters
	config = ConfigParser.SafeConfigParser()
	config.optionxform = str
	config.read(fn)
	def getVal(val, defaultVal):
		if isinstance(defaultVal, int):
			return int(val)
		elif isinstance(defaultVal, float):
			return float(val)
		elif isinstance(defaultVal, str):
			return str(val)
		else:
			errAbort("Unsure how to convert %s" % val)

	for section in AllFilters.keys():
		if config.has_section(section):
			for var, val in config.items(section):
				if var not in AllFilters[section]:
					errAbort("%s filter %s is not a valid filter." % (section, var))
				AllFilters[section][var] = getVal(val, AllFilters[section][var])
			
	return config

def main():
	''' Main function for Filter '''
	global KEEP_TMP_FILES, AllFilters
	args = parseArgv()
	mutations = []
	configFn, pointMutFn, indelFn, outFn = args

	config = configure(configFn)
	filterPointMutations(pointMutFn, mutations)
	print "Point mutations: %d" % len(mutations)
	filterIndels(indelFn, mutations)
	print "Total mutations: %d" % len(mutations)
	annoVarInputFn = writePreAnnovarMuts(mutations)
	annoVarResults = getAnnoVarResults(annoVarInputFn, str(config.get('General', 'genomeBuild')), str(config.get('AnnoVarInputs', 'executable')), str(config.get('AnnoVarInputs', 'dbDir')), str(config.get('AnnoVarInputs', '1kgDBs')), str(config.get('AnnoVarInputs', 'snpDBs')), AllFilters['AnnoVarAnnotationFilters'])
	print "Annovar mutations: %d" % len(annoVarResults)
	writeMutations(outFn, mutations, annoVarResults)

	# Clean up
	if not KEEP_TMP_FILES:
		for tmpFn in glob.glob(annoVarInputFn + "*"):
			os.remove(tmpFn)

def parseArgv():
	''' Parse the command line options and return the arguments '''
	global NUM_ARGS, KEEP_TMP_FILES, TMP_DIR
	try:
		opts, args = getopt.gnu_getopt(sys.argv[1:], "h?", ["help","keepTmpFiles","tmp="])
	except getopt.error, msg:
		print msg
		print __doc__
		sys.exit(-1)

	if len(args) != NUM_ARGS:
		print "Filter.py requires %d arguments, %d given." % (NUM_ARGS, len(args))
		print __doc__
		sys.exit(-1)

	for o,a in opts:
		if o in ("-h", "-?", "--help"):
			print __doc__
			sys.exit(0)
		elif o == "--keepTmpFiles":
			KEEP_TMP_FILES = True
		elif o == "--tmp":
			TMP_DIR = a
		else:
			print "Invalid option:", o
			print __doc__
			sys.exit(-1)

	return args


def filterPointMutations(mutFn, retList):
	''' Append valid point mutations to the retList '''
	global AllFilters
	try:
		f = open(mutFn)
	except IOError:
		errAbort("Point mutation file %s does not exist." % mutFn)

	version = f.readline().strip()
	rawHeader = f.readline().strip()
	numCols = validateMutectorFile(version, rawHeader.replace('\t', ' '))

	snpRemovalFilters = AllFilters['SNPRemovalFilters']
	snpAnnotationFilters = AllFilters['SNPAnnotationFilters']
	print "Screening SNP mutations with the following parameters:"
	print "  ", "For removal of the call:"
	for k in sorted(snpRemovalFilters.keys()):
		print "    ", k, ":", snpRemovalFilters[k]
	print "  ", "For annotation of the call:"
	for k in sorted(snpAnnotationFilters.keys()):
		print "    ", k, ":", snpAnnotationFilters[k]

	# Massage removal filters to be more efficient
	snpRemovalFilters['muTectJudgmentsAllowed'] = set(snpRemovalFilters['muTectJudgmentsAllowed'].split(','))

	# Loop through each successive line in the file
	while True:
		line = f.readline()
		if not line:
			break
		mutation = MuTectorLine(line, snpAnnotationFilters)
		if mutation.shouldBeKept(snpRemovalFilters):
			retList.append(mutation)


def filterIndels(indelFn, retList):
	''' Append valid indel mutations to the retList '''
	global AllFilters
	try:
		f = open(indelFn)
	except IOError:
		errAbort("Indel mutation file %s does not exist." % indelFn)

	inDataLine = False
	tumorName = ""; normalName = ""

	indelRemovalFilters = AllFilters['IndelRemovalFilters']
	indelAnnotationFilters = AllFilters['IndelAnnotationFilters']
	print "Screening indels with the following parameters:"
	print "  ", "For removal of the call:"
	for k in sorted(indelRemovalFilters.keys()):
		print "    ", k, ":", indelRemovalFilters[k]
	print "  ", "For annotation of the call:"
	for k in sorted(indelAnnotationFilters.keys()):
		print "    ", k, ":", indelAnnotationFilters[k]

	# Loop through each successive line in the file
	while True:
		line = f.readline()
		if not line:
			break
		if line.startswith("##"):
			if inDataLine:
				errAbort("Invalid format for .vcf, should not have any comments after initial headers: %s" % line)
		elif line.startswith("#"):
			if inDataLine:
				errAbort("Invalid format for .vcf, should not have any comments after initial headers: %s" % line)
			spaceLine = line.strip().replace('\t', ' ')
			if not spaceLine.startswith(SomaticIndelRequiredColumns):
				errAbort("Expecting .vcf to have different columns than: %s" % line)
			tumorName, normalName = line.strip().split('\t')[-2:]
		else:
			inDataLine = True
			mutation = SomaticIndelLine(line, tumorName, normalName, indelAnnotationFilters)
			if mutation.shouldBeKept(indelRemovalFilters):
				retList.append(mutation)

def writePreAnnovarMuts(mutations):
	''' Write out the mutations that pass all filters in an AnnoVar-acceptable format '''
	global TMP_DIR
	currDateTimeStr = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
	prefix = "AnnoVar.%s." % currDateTimeStr
	ntfNum, fn = tempfile.mkstemp(text=True, dir=TMP_DIR, prefix=prefix)
	f = open(fn, "w")
	for mutation in mutations:
		f.write("%s\n" % mutation.AnnoVarStr())
	f.close()
	return fn


def writeMutations(outFn, retval, annovars):
	''' Write out all valid mutations '''
	global AllFilters
	assert len(retval) == len(annovars)
	f = open(outFn, "w")
	f.write("#gene\tcontig\tposition\tref_allele\talt_allele\tnucleotide\tprotein\tcontext\ttype\ttumor_name\tnormal_name\tscore\tpower\ttumor_power\tnormal_power\ttotal_pairs\timproper_pairs\tmap_Q0_reads\tcontaminant_fraction\tcontaminant_lod\tt_ref_count\tt_alt_count\tn_ref_count\tn_alt_count\ttumor_variant_freq\tnormal_variant_freq\taccession\texon\tknown_variant_status\talgorithm\tourJudgment\tourReason\n")
	
	def getOurFinalJudgment(j1, j2):
		if j1 == "yes" and j2 == "yes":
			return "yes"
		return "no"

	def joinText(r1, r2, defaultText="", joinerText=","):
		if r1 == defaultText:
			return r2
		if r2 == defaultText:
			return r1
		return r1 + joinerText + r2

	annoVarRemovalFilters = AllFilters['AnnoVarRemovalFilters']

	for idx, m in enumerate(retval):
		# We've already filtered the MuTect and SomaticIndel data based on reads
		# and variant frequency.  Now do the final filters to make sure it's a type
		# of mutation we care about reporting.
		a = annovars[idx]

		# Check that everything looks okay
		assert (a.chrom == m.contig) and (a.startPos == m.position) and (a.ref_allele == m.ref_allele) and (a.alt_allele == m.alt_allele)

		if a.shouldBeKept(annoVarRemovalFilters): # and getOurFinalJudgment(m.ourJudgment,a.ourJudgment) == "yes":
			if m.algo == SomaticIndelAlgorithm:
				f.write("\t".join([a.genes, a.chrom, str(a.startPos), a.ref_allele, a.alt_allele, a.nucleotide, a.protein, a.context, a.exonMutTypeStr(), m.tumor_name, m.normal_name, str(m.score), str(m.power), str(m.tumor_power), str(m.normal_power), str(m.t_total_depth()+m.n_total_depth()), str(m.improper_pairs), str(m.map_Q0_reads()), "0", NA, str(m.t_ref_count()), str(m.t_alt_count()), str(m.n_ref_count()), str(m.n_alt_count()), str(m.t_var_freq()), str(m.n_var_freq()), a.accession, a.exon, joinText(m.getKnownVariantStatus(), a.getKnownVariantStatus(), defaultText="NOVEL", joinerText="_"), m.algo, getOurFinalJudgment(m.ourJudgment, a.ourJudgment), joinText(m.ourJudgmentReasons, a.ourJudgmentReasons)]) + '\n')
			else:
				f.write("\t".join([a.genes, a.chrom, str(a.startPos), a.ref_allele, a.alt_allele, a.nucleotide, a.protein, a.context, a.exonMutTypeStr(), m.tumor_name, m.normal_name, str(m.score), str(m.power), str(m.tumor_power), str(m.normal_power), str(m.total_pairs), str(m.improper_pairs), str(m.map_Q0_reads), str(m.contaminant_fraction), str(m.contaminant_lod), str(m.t_ref_count), str(m.t_alt_count), str(m.n_ref_count), str(m.n_alt_count), str(m.t_var_freq()), str(m.n_var_freq()), a.accession, a.exon, joinText(m.getKnownVariantStatus(), a.getKnownVariantStatus(), defaultText="NOVEL", joinerText="_"), m.algo, getOurFinalJudgment(m.ourJudgment, a.ourJudgment), joinText(m.ourJudgmentReasons, a.ourJudgmentReasons)]) + '\n')
	f.close()


if __name__ == '__main__':
	sys.exit(main())
