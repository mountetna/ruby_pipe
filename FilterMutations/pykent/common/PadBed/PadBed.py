#!/usr/bin/env python

"""
PadBed.py - A class to place padding around Bed elements.

Usage: PadBed.py chromSizes input.bed paddedOutput.bed
  where
  chromSizes        is either a .2bit file of the genome used or an explicit chrom.sizes file
  input.bed         is the input file that should be padded
  paddedOutput.bed  is the output file with padding added

Options:
  --twoBit        chromSizes file is a .2bit file
  --sizes         chromSizes file is a .sizes file
  --padding=N     Pad both upstream and downstream by N base pairs (default=0)
  --upstream=N    For stranded Beds, pad upstream by N base pairs (default=0)
  --downstream=N  For stranded Beds, pad downstream by N base pairs (default=0)
"""

import sys
import os
import getopt
from pykent.common.Sanity import errAbort
from pykent.common.GenomeUtil import getChromSizesDict, getChromSizesDictFromText
import pykent.common.Bed as Bed

NUM_ARGS = 3
IS_TWO_BIT = False
IS_SIZES = False
PADDING = -1
UPSTREAM = -1
DOWNSTREAM = -1

def main():
	''' Main function for PadBed '''
	global IS_TWO_BIT, IS_SIZES
	args = parseArgv()
	chromSizesFn, inFn, outFn = args

	# Check if chromSizesFn type has been set explicitly by arguments, if not, try to guess
	if (not IS_TWO_BIT) and (not IS_SIZES):
		if chromSizesFn.endswith(".2bit"):
			IS_TWO_BIT = True
		elif chromSizesFn.endswith(".sizes"):
			IS_SIZES = True
		else:
			errAbort("Unknown file type for %s.  Please either end in '.2bit' or '.sizes' or set the appropriate flag." % chromSizes)

	if IS_TWO_BIT:
		chromSizesDict = getChromSizesDict(chromSizesFn)
	else:
		chromSizesDict = getChromSizesDictFromText(chromSizesFn)

	padBed(chromSizesDict, inFn, outFn)

def parseArgv():
	''' Parse the command line options and return the arguments '''
	global NUM_ARGS, IS_TWO_BIT, IS_SIZES, PADDING, UPSTREAM, DOWNSTREAM
	try:
		opts, args = getopt.gnu_getopt(sys.argv[1:], "h?", ["help", "twoBit", "sizes", "padding=", "upstream=", "downstream="])
	except getopt.error, msg:
		print msg
		print __doc__
		sys.exit(-1)

	if len(args) != NUM_ARGS:
		print "PadBed.py requires %d arguments, %d given." % (NUM_ARGS, len(args))
		print __doc__
		sys.exit(-1)

	for o,a in opts:
		if o in ("-h", "-?", "--help"):
			print __doc__
			sys.exit(0)
		elif o == "--twoBit":
			IS_TWO_BIT = True
		elif o == "--sizes":
			IS_SIZES = True
		elif o == "--padding":
			PADDING = int(a)
			if PADDING < 0:
				errAbort("Padding must be a non-negative integer: %d" % PADDING)
		elif o == "--upstream":
			UPSTREAM = int(a)
			if UPSTREAM < 0:
				errAbort("Upstream padding must be a non-negative integer: %d" % UPSTREAM)
		elif o == "--downstream":
			DOWNSTREAM = int(a)
			if DOWNSTREAM < 0:
				errAbort("Downstream padding must be a non-negative integer: %d" % DOWNSTREAM)

	if IS_TWO_BIT and IS_SIZES:
		errAbort("Must specify only one of '--twoBit', '--sizes' flags.")

	if PADDING >= 0:
		if (UPSTREAM >= 0) or (DOWNSTREAM >= 0):
			errAbort("Can specify either global padding (using '--padding') or strand-specific padding (using '--up/downstream'), but not both.")

	return args

def padBed(sizes, inFn, outFn):
	''' Perform the actual padding of the Beds. '''
	global PADDING, UPSTREAM, DOWNSTREAM
	isStrandedPadding = (PADDING < 0)
	if not isStrandedPadding:
		UPSTREAM = PADDING
		DOWNSTREAM = PADDING
	else:
		UPSTREAM = max(0, UPSTREAM)
		DOWNSTREAM = max(0, DOWNSTREAM)

	f = open(inFn)
	g = open(outFn, "w")
	validStrands = ('+', '-')
	for line in f:
		b = Bed.Bed(line, chromSizes=sizes)
		if isStrandedPadding:
			if (not hasattr(b, 'strand')) or (b.strand not in validStrands):
				errAbort("Strand-specific padding requires all input elements have valid strands: %s" % b)
			elif b.strand == '+':
				b.chromStart = max(0, b.chromStart-UPSTREAM)
				b.chromEnd = min(sizes[b.chrom], b.chromEnd+DOWNSTREAM)
			elif b.strand == '-':
				b.chromStart = max(0, b.chromStart-DOWNSTREAM)
				b.chromEnd = min(sizes[b.chrom], b.chromEnd+UPSTREAM)
			else:
				errAbort("Should not happen, programmer error.")
		else:	
			b.chromStart = max(0, b.chromStart-UPSTREAM)
			b.chromEnd = min(sizes[b.chrom], b.chromEnd+DOWNSTREAM)

		g.write("%s\n" % b)
	f.close()
	g.close()

if __name__ == '__main__':
	sys.exit(main())

