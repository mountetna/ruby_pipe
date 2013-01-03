#!/usr/bin/env python

"""
A utility class to provide methods operating on Bed entries
"""

from operator import attrgetter

import Bed
from Sanity import errAbort

def sortByChromStartEnd(lst):
	''' Sort a list of Bed elements in place by chromosome, start, and end position '''
	lst[:] = sorted(lst, key=attrgetter('chrom', 'chromStart', 'chromEnd'))

def bedReadFromFile(fn, chromSizes=None):
	''' Read a file of Bed regions into a list of Beds and return it '''
	f = open(fn)
	retval = [Bed.Bed(line, chromSizes=chromSizes) for line in f.readlines()]
	f.close()
	return retval

def bedChromDictFromFile(fn, chromSizes=None):
	''' Return a dictionary keyed by chromosome with values being all Bed regions on that chrom (sorted) '''
	retval = {}
	bedList = bedReadFromFile(fn, chromSizes)
	for b in bedList:
		retval[b.chrom] = retval.get(b.chrom, []) + [b]
	for v in retval.values():
		sortByChromStartEnd(v)
	return retval

def isNonContiguous(lst, isSorted=True):
	''' Return True iff Beds are not overlapping or contiguous.  Modifies list if sorted is not True '''
	if not isSorted:
		sortByChromStartEnd(lst)

	lstLen = len(lst)
	if lstLen <= 1:
		return True

	for i in xrange(1, lstLen):
		prev = lst[i-1]
		curr = lst[i]
		if curr.chrom < prev.chrom:
			errAbort("Bed list is not sorted.")

		elif curr.chrom == prev.chrom:
			if (curr.chromStart < prev.chromStart) or ((curr.chromStart == prev.chromStart) and (curr.chromEnd < prev.chromEnd)):
				errAbort("Bed list is not sorted.")

			if curr.chromStart <= prev.chromEnd:
				print prev.chrom, prev.chromStart, prev.chromEnd
				return False
	return True

def isNonContiguousDict(d, isSorted=True):
	''' Return True iff all Beds in the dictionary are entirely non-contiguous.  Modifies unsorted lists if required. '''
	retval = True
	for bedList in d.values():
		if not isNonContiguous(bedList, isSorted):
			retval = False
	return retval

def isBedListSorted(lst):
	''' Return True iff lst is sorted by chrom, chromStart, chromEnd '''
	lstLen = len(lst)
	for i in xrange(1, lstLen):
		prev = lst[i-1]
		curr = lst[i]
		if curr.chrom < prev.chrom:
			return False
		elif curr.chrom == prev.chrom:
			if (curr.chromStart < prev.chromStart) or ((curr.chromStart == prev.chromStart) and (curr.chromEnd < prev.chromEnd)):
				return False
	return True

def isNonOverlapping(lst, isSorted=True):
	''' Return True iff Beds are entirely non-overlapping.  Modifies list if sorted is not True '''
	if not isSorted:
		sortByChromStartEnd(lst)

	lstLen = len(lst)
	if lstLen <= 1:
		return True

	for i in xrange(1, lstLen):
		prev = lst[i-1]
		curr = lst[i]
		if curr.chrom < prev.chrom:
			errAbort("Bed list is not sorted.")

		elif curr.chrom == prev.chrom:
			if (curr.chromStart < prev.chromStart) or ((curr.chromStart == prev.chromStart) and (curr.chromEnd < prev.chromEnd)):
				errAbort("Bed list is not sorted.")

			if curr.chromStart < prev.chromEnd:
				return False
	return True

def isNonOverlappingDict(d, isSorted=True):
	''' Return True iff all Beds in the dictionary are entirely non-overlapping. Modifies unsorted lists if required. '''
	retval = True
	for bedList in d.values():
		if not isNonOverlapping(bedList, isSorted):
			retval = False
	return retval

def isEncompassedByBed(chrom, start, end, bedList, bedListSorted=False):
	''' Return True iff the chrom, start, end is encompassed by a single element within bedList. '''
	if not bedListSorted:
		sortByChromStartEnd(bedList)

	isEncompassed = False
	for bed in bedList:
		if chrom < bed.chrom:
			break
		elif chrom > bed.chrom:
			continue
		else:
			if start < bed.chromStart:
				break
			elif start > bed.chromEnd:
				continue
			else:
				if end <= bed.chromEnd:
					isEncompassed = True
					break
	return isEncompassed

def getOverlappingRegionDict(bd1, bd2, debug=True):
	''' Return a dictionary of overlapping regions within two dictionaries of Bed regions '''
	retval = {}
	for chrom, bedList1 in bd1.items():
		bedList2 = bd2.get(chrom, [])

		if bedList2 == []:
			continue

		if debug:
			if not isNonContiguous(bedList1):
				errAbort("Calculating overlapping regions must have non-contiguous input elements.")
			if not isNonContiguous(bedList2):
				errAbort("Calculating overlapping regions must have non-contiguous input elements.")

		b1Len = len(bedList1); b1Idx = 0
		b2Len = len(bedList2); b2Idx = 0

		while (b1Idx < b1Len) and (b2Idx < b2Len):
			b1Curr = bedList1[b1Idx]
			b2Curr = bedList2[b2Idx]
			assert b1Curr.chrom == b2Curr.chrom
			maxStart = max(b1Curr.chromStart, b2Curr.chromStart)
			minEnd = min(b1Curr.chromEnd, b2Curr.chromEnd)
			if maxStart < minEnd:
				retval[b1Curr.chrom] = retval.get(b1Curr.chrom,[]) + [Bed.Bed("%s\t%d\t%d" % (b1Curr.chrom,maxStart,minEnd))]

			if b1Curr.chromEnd < b2Curr.chromEnd:
				b1Idx += 1
			elif b1Curr.chromEnd > b2Curr.chromEnd:
				b2Idx += 1
			else:
				b1Idx += 1; b2Idx += 1
	return retval

