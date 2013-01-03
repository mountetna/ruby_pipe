#!/usr/bin/env python

"""
File of sanity checking functions and operations for use by all source code.
"""

import sys

def errAbort(msg, exitCode=255):
	''' Print error message and exit entirely '''
	print "Err:", msg
	sys.exit(exitCode)

def canBeNum(n):
	''' Return True iff input variable is an int, float, or string that represents one '''
	if isInt(n) or isFloat(n):
		return True
	try:
		stripped = str(float(n))
		return True
	except ValueError:
		return False

def canBeInt(n):
	''' Return True iff input variable is an int or a string representing one '''
	if isInt(n):
		return True
	try:
		stripped = str(int(n))
		return True
	except ValueError:
		return False

def isInt(n):
	''' Return True iff input variable is an integer '''
	return isinstance(n, int)

def isNonNegInt(n):
	''' Return True iff input variable is a non-negative integer '''
	return isInt(n) and (n >= 0)

def isNonPosInt(n):
	''' Return True iff input variable is a non-positive integer '''
	return isInt(n) and (n <= 0)

def isPosInt(n):
	''' Return True iff input variable is a positive integer '''
	return isInt(n) and (n > 0)

def isNegInt(n):
	''' Return True iff input variable is a negative integer '''
	return isInt(n) and (n < 0)

def isFloat(n):
	''' Return True iff input variable is a float '''
	return isinstance(n, float)

def isNonNegFloat(n):
	''' Return True iff input variable is a non-negative float '''
	return isFloat(n) and (n >= 0)

def isNonPosFloat(n):
	''' Return True iff input variable is a non-positive float '''
	return isFloat(n) and (n <= 0)

def isPosFloat(n):
	''' Return True iff input variable is a positive float '''
	return isFloat(n) and (n > 0)

def isNegFloat(n):
	''' Return True iff input variable is a negative float '''
	return isFloat(n) and (n < 0)

def compare(a, b):
	if (type(a) != type(b)) and not (isinstance(a, (int, float, long)) and isinstance(b, (int, float, long))):
		errAbort("Type mismatch between objects to compare: %s (%s), %s (%s)" % (str(a), str(type(a)), str(b), str(type(b))))
	if a < b:
		return -1
	elif a > b:
		return 1
	return 0

def inRange(n, minVal, maxVal):
	''' Return True iff n is within the range of [minVal, maxVal] '''
	if ((type(minVal) == type(maxVal)) and (type(minVal) == type(n))) or \
		(isinstance(minVal, (int, float, long)) and isinstance(maxVal, (int, float, long)) and isinstance(n, (int, float, long))):
		return (minVal <= n) and (n <= maxVal)
	else:
		errAbort("minVal and maxVal are not of the same type: %s, %s" % (str(minVal), str(maxVal)))

def roughlyEqual(f1, f2, epsilon=1e-6):
	''' Return True iff floats f1 and f2 are within epsilon of each other '''
	return (abs(f1-f2) < epsilon)

