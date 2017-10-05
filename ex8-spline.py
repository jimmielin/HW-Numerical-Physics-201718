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
# Now Playing: Remember the Name - Fort Minor
#
# Reticulating Splines requires the solving of a tri-diagonal matrix
# using whatever method of our choice. I have chosen Thomas' Algorithm,
# which is similar to 追赶法 because its easier and more lightweight
# to implement.
# The implementation here is similar to a previous
# masterpiece of mine written in Haskell,
# https://github.com/jimmielin/HW-Numerical-Analysis-Spring-2016/
#  blob/master/20160305-extrapolation/splines.hs
# but re-implemented in Python.
####################################################################

import math
import time

# Include my own lightweight matrix library & its prototype
from pylitematrix.pylitematrix import Matrix, mp

####################################################################
# The Tridiagonal Matrix Algorithm "Thomas Algorithm"
# https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
#
# An algorithm stable for solving tri-diagonal, positive semi-definite
# or diagonally dominant matrices
# in O(n) instead of O(n^2) as Gauss.

# solveTri(M, b) = x solves Mx=b
# when M is tri-diagonal, positive semi-definite || diagonal dominant
#
# Notes of the author (10/05/2017):
# 1. This implementation IS space-optimal as the matrix is decomposed in
# the first step and stored in an internal representation similar to the
# course description (but since we are not symmetric, we use arrays)
# 
def solveTri(M, b):
	# Extract the {c_i}, {b_i} and {a_i} in the NxN matrix
	n = M.nrows()
	if(M.ncols() != n):
		raise ValueError("Tridiagonal Matrix to-solve must be square.")

	# The arrays below have been padded so the indeces are mathematically sound.
	# One doesn't have to do this, but it makes for more coding convenience.
	A = [0, 0] # lower diagonal A, i=2,3..n
	B = [0] # diagonal B, i=1,2,..n
	C = [0] # upper diagonal C, i=1,2,..n-1
	D = b.toList() # right-b vector d1, d2, ... n
	D.insert(0, 0) # left-padd, * minor performance concern

	# The solution array is *NOT* padded, as it'll be used for
	# fromList_ by the Matrix prototype
	# It's even reversed as part of back-substitution
	# so watch out for the index.
	X = [None] * n

	# .elem is now O(1) so you can use that for safety,
	# but I'll tap into .internal represenation for speed and less overhead
	# since we already checked the indexes above
	for i in range(1, n+1):
		if(i != 1): # make sure i=1 case does not run for A
			A.append(M.internal[i-1][i-2])
		if(i != n): # make sure i=n case does not run for C
			C.append(M.internal[i-1][i])
		B.append(M.internal[i-1][i-1])

	# Now perform the Thomas Algorithm Forward Sweep
	C[1] = C[1] / B[1]
	D[1] = D[1] / B[1]
	for i in range(2, n+1):
		if(i != n): # make sure i=n case does not run for C=..n-1
			C[i] = C[i] / (B[i] - A[i] * C[i-1])
		D[i] = (D[i] - A[i] * D[i-1]) / (B[i] - A[i] * C[i-1])

	# Now do a back-subsitution...
	X[0] = D[n]
	for i in range(2, n+1):
		j = n+1-i # back indexes
		print(j)
		X[i-1] = D[j] - C[j] * X[i-2]

	X.reverse()
	print(X)

# Some testing data.
#testtrim = mp.fromList_(4, 4, [1,5,0,0,8,2,6,0,0,9,3,7,0,0,10,4])
#testtrib = mp.fromList_(4, 1, [1,1,1,1])
#print(testtrim)
