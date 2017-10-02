####################################################################
# Computational Physics, 2017-18 Sem1
# HW-1 Ex-8
#
# (c) 2017 Haipeng Lin <linhaipeng@pku.edu.cn>
# All Rights Reserved.
#
# This program is written as Homework for the Computational Physics
# Course in Peking University, School of Physics.
# NO WARRANTIES, EXPRESS OR IMPLIED, ARE OFFERED FOR THIS PRODUCT.
# Copying is strictly prohibited and its usage is solely restricted
# to uses permitted by the author and homework grading.
#
# This program is Python 3 (3.6.2 on darwin)
# compatible, and does not support Python 2.
#
# Now Playing: Forgettable - Project 46
####################################################################

import math
import time

# Include my own lightweight matrix library
from pylitematrix.pylitematrix import Matrix

# Construct the multi-diagonal matrix of Ex. 7 by definition,
# directly using an internal multi-diagonal matrix representation C
# (formula No. 67, Chapter 2)
# that is appropriate in space-complexity.
#
# every row of C is a diagonal of A, with the top-left zero-padded.
# a(i, j) translates to c(i, j - i + m + 1)
# (when i<=m+1, j=1,2,3...i) (when i>m+1, j=i-m,...,i)
# where A is a symmetric, positive-definite, multi-diagonal matrix
# with a band-width of 2m+1 (def. as |i-j|>m => a(ij)=0), m is half width
#
# then C is (spatially) a (n)x(m+1) matrix
#
# A is N-deg, m=2
def reprMatrC(n):
	# another rant:
	# the reason this lambda is so messy is because Guido van Rossum decided
	# that multi-line lambdas would be a "Rube Goldberg Contraption" and
	# all proposed solutions would be non-pythonic.
	# whatever that means.
	# one has to have a indentation-based, lambda-lightweight, side-effect
	# "bad" programming language for the masses in this world, right?
	# Python is in some ways even worse in design than PHP.

	# a(1,1)=a(n,n)=5 => c(1,3)=c(n,3)=5
	# a(i,i)=6 otherwise => c(i,3=6)
	# a(i,i-1)=4 => c(i,2)=4
	# a(i,i-2)=1 => c(i,1)=1
	return Matrix(n, 3, lambda i, j: 1 if j == 1 else (4 if j == 2 else (5 if (i == 1 or i == n) else 6)))

def convertAtoCcoords(i, j): # unsafe, for m = 2 only
	return (i, j-i+3)

def convertCtoAcoords(i, j): # unsafe, for m = 2 only
	return (i, j+i-3)

# using the sqrt-method calculate matrix A (repr C)'s Cholesky decomposition
def sqrtMethodL(C):
	