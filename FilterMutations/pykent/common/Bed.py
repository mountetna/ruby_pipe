#!/usr/bin/env python

from Sanity import *

class Bed:

	MIN_FIELDS = 3
	MAX_FIELDS = 15
	STRAND_FIELDS = ('+', '-', '.')
	
	def __init__(self, line, chromSizes=None, sep='\t', intScore=True):
		''' Create a Bed structure from a line of text '''
		if line.endswith('\n'):
			line = line[:-1]
		splitLine = line.split(sep)
		numFields = len(splitLine)

		if not inRange(numFields, Bed.MIN_FIELDS, Bed.MAX_FIELDS):
			errAbort("Bed line must have between %d and %d fields, not %d" % (Bed.MIN_FIELDS, Bed.MAX_FIELDS, numFields))

		self.numFields = numFields

		chrom, chromStart, chromEnd = splitLine[0:3]
		self.chrom = chrom
		self.chromStart = int(chromStart)
		self.chromEnd = int(chromEnd)
		if self.chromStart < 0:
			errAbort("chromStart (%d) must be non-negative." % self.chromStart)
		if self.chromStart > self.chromEnd:
			errAbort("chromStart (%d) must be <= chromEnd (%d)" % (self.chromStart, self.chromEnd))

		if chromSizes != None:
			if self.chrom not in chromSizes:
				errAbort("chrom (%s) not present in chromSizes dictionary." % self.chrom)
			if self.chromEnd > chromSizes[self.chrom]:
				errAbort("chromEnd (%d) larger than the chromSizes listing for this chromosome (%d)" % (self.chromEnd, chromSizes[self.chrom]))

		if self.numFields > 3:
			self.name = splitLine[3]

			if self.numFields > 4:
				self.score = splitLine[4]
				if not (canBeInt(self.score) or (canBeNum(self.score) and not intScore)):
					req = " number"
					if intScore:
						req = "n integer"
					errAbort("Score field must be a%s: %s" % (req, self.score))
				if canBeInt(self.score):
					self.score = int(self.score)
				else:
					self.score = float(self.score)

				if self.numFields > 5:
					self.strand = splitLine[5]
					if self.strand not in Bed.STRAND_FIELDS:
						errAbort("Strand must be one of the supported types: %s" % ",".join(Bed.STRAND_FIELDS))

					if self.numFields > 6:
						if self.numFields < 8:
							errAbort("If thickStart is specified, thickEnd must also be specified.")
						self.thickStart = int(splitLine[6])
						self.thickEnd = int(splitLine[7])
						if self.thickStart < self.chromStart:
							errAbort("thickStart (%d) must be >= chromStart (%d)" % (self.thickStart, self.chromStart))
						if self.thickStart > self.thickEnd:
							errAbort("thickStart (%d) must be <= thickEnd (%d)" % (self.thickStart, self.thickEnd))
						if self.thickEnd > self.chromEnd:
							errAbort("thickEnd (%d) must be <= chromEnd (%d)" % (self.thickEnd, self.chromEnd))

						if self.numFields > 8:
							self.itemRgb = splitLine[8]

							if self.numFields > 9:
								if self.numFields < 12:
									errAbort("If blockCount is specified, blockSizes and blockStarts must also be specified")
								self.blockCount = int(splitLine[9])
								self.blockSizes = [int(x) for x in splitLine[10].rstrip(',').split(',')]
								self.blockStarts = [int(x) for x in splitLine[11].rstrip(',').split(',')]

								if (self.blockCount != len(self.blockSizes)) or (self.blockCount != len(self.blockStarts)):
									errAbort("blockCount does not match the written number of blockSizes or blockStarts")

								if len(filter(isNonPosInt, self.blockSizes)) > 0:
									errAbort("Block sizes must be positive integers %s" % self.blockSizes)

								if len(filter(isNegInt, self.blockStarts)) > 0:
									errAbort("Block starts must be non-negative integers %s" % self.blockStarts)
								
								for n in xrange(1, self.blockCount):
									if self.blockStarts[n-1]+self.blockSizes[n-1] > self.blockStarts[n]:
										errAbort("Overlapping blocks not allowed within Bed: %s\t%s" % (",".join(map(str, self.blockSizes)), ",".join(map(str, self.blockStarts))))

								if self.blockStarts[-1]+self.blockSizes[-1] > self.chromEnd:
									errAbort("Blocks must not extend beyond the end of the Bed %d, %d+%d" % (self.chromEnd, self.blockSizes[-1], self.blockStarts[-1]))

								if self.numFields > 12:
									self.rest = sep.join(splitLine[12:])

	def __str__(self):
		''' Return a string representation of the Bed element '''
		retval = "%s\t%d\t%d" % (self.chrom, self.chromStart, self.chromEnd)
		if self.numFields > 3:
			retval += "\t%s" % self.name
			if self.numFields > 4:
				if isInt(self.score):
					retval += "\t%d" % self.score
				else:
					retval += "\t%f" % self.score
				if self.numFields > 5:
					retval += "\t%s" % self.strand
					if self.numFields > 6:
						retval += "\t%d\t%d" % (self.thickStart, self.thickEnd)
						if self.numFields > 8:
							retval += "\t%s" % self.itemRgb
							if self.numFields > 9:
								retval += "\t%d\t%s\t%s" % (self.blockCount, ",".join(map(str, self.blockSizes))+",", ",".join(map(str, self.blockStarts))+",")
								if self.numFields > 12:
									retval += "\t%s" % self.rest
		return retval

	def __cmpByChromStartEnd__(self, other):
		''' Compare Bed elements based on chrom, chromStart, and chromEnd values '''
		if self.chrom < other.chrom:
			return -1
		elif self.chrom > other.chrom:
			return 1

		if self.chromStart < other.chromStart:
			return -1
		elif self.chromStart > other.chromStart:
			return 1

		if self.chromEnd < other.chromEnd:
			return -1
		elif self.chromEnd > other.chromEnd:
			return 1

		return 0

