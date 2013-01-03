#!/usr/bin/env python

"""
A utility class for functions on genomes
"""

import subprocess
from Sanity import errAbort

def getChromSizesDict(fn):
	''' Return a dictionary keyed by chromosome with value as chrom size based on input .2bit file '''
	retval = {}
	cmd = "twoBitInfo %s stdout" % fn
	proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
	output = str(proc.communicate()[0])
	for line in output.strip().split('\n'):
		chrom, chromSize = line.strip().split('\t')
		retval[chrom] = int(chromSize)
	return retval

def getChromSizesDictFromText(textFn):
	''' Return a dictionary keyed by chromosome with value as chrom size based on input chrom.sizes file '''
	retval = {}
	f = open(textFn)
	for line in f:
		chrom, chromSize = line.strip().split('\t')
		retval[chrom] = int(chromSize)
	return retval

def getDNA(chrom, start, end, fn, noMask=False):
	''' Return the DNA associated with the BED-style position as a single long string '''
	maskFlag = ""
	if noMask:
		maskFlag = " -noMask"
	cmd = ("twoBitToFa %s stdout" % fn) + maskFlag + " -seq="
	proc = subprocess.Popen(cmd + "%s:%d-%d" % (chrom, start, end), shell=True, stdout=subprocess.PIPE)
	dna = str(proc.communicate()[0]).strip()
	dna = dna.split('\n')
	if len(dna) < 2:
		errAbort("Must be at least a header line and one line of DNA.")
	tmp = "".join(dna[1:])
	return "".join(tmp.split())

