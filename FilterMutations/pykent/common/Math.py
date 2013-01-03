#!/usr/bin/env python

"""
A utility file for many mathematical functions
"""
from math import sqrt

def mean(lst):
	''' Return the arithmetic mean of a list of numbers, or None if the list is empty '''
	if len(lst) == 0:
		return None
	return sum(lst)*1./len(lst)

def stdev(lst):
	''' Return the sample standard deviation of a list of numbers, or None if the list is empty '''
	numElts = len(lst)
	if numElts == 0:
		return None
	if numElts == 1:
		return 0

	m = mean(lst)
	squaredErrors = [(x-m)**2 for x in lst]
	sampSum = sum(squaredErrors)*1./(numElts-1)
	return sqrt(sampSum)

def corr(x,y):
	''' Return the correlation coefficient between the data sets '''
	numElts = len(x)
	assert len(y) == numElts
	if numElts == 0:
		return None

	xMean = mean(x); yMean = mean(y)
	xStd = stdev(x); yStd = stdev(y)

	numer = sum([(x[i]-xMean)*(y[i]-yMean) for i in xrange(numElts)])
	return numer*1./((numElts-1)*xStd*yStd)

