# $ RELEASE $
# $ 201711060211Z $ rel01
# $ Signed-Off-By: Haipeng Lin <jimmie.lin@gmail.com>
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

from functools import partial

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
# 2. Watch out for the indexes. Some have been padded for mathematical
# understanding, but some have not for implementation reasons.
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
        X[i-1] = D[j] - C[j] * X[i-2]

    X.reverse()

    return mp.fromList_(n, 1, X)

# Some testing data.
#testtrim = mp.fromList_(4, 4, [1,5,0,0,8,2,6,0,0,9,3,7,0,0,10,4])
#testtrib = mp.fromList_(4, 1, [1,1,1,1])
#print(solveTri(testtrim, testtrib))

####################################################################
# Reticulator of Splines
# “三次样条差值是一种分段三次多项式插值, 在每个插值节点处比分段三次Hermite更光滑,
#  具有二阶连续导数, 而且不需要节点导数信息.”
# 据王鸣老师所述，在早期工程学上三次样条插值的实现方法就是在坐标纸上画点, 然后用
# 一个铁条固定在节点处弯曲，从而获得自然边界条件。
#
# 王教授于2016年夏因病逝世, 我们成为最后一批《计算方法B》的学生.
# R.I.P.
#

# Base Hermite Functions alpha, beta
# "分段三次Hermite插值的基函数"
#
# for partial applications we may need functools->partial, but this implementation
# is so much lamer than Haskell...
aL = lambda a, b, x: (1 + 2*(x-a)/(b-a)) * ((x-b)/(a-b))**2
aR = lambda a, b, x: (1 + 2*(x-b)/(a-b)) * ((x-a)/(a-b))**2
bL = lambda a, b, x: (x-a)*((x-b)/(a-b))**2
bR = lambda a, b, x: (x-b)*((x-a)/(a-b))**2

# Hermite Functions alpha_i, beta_i
# Generated Functions (Partially Applied Lambdas) that use aLRbLR to
# generate the respective alpha/beta functions for extrapolation
# According to 《计算方法》Page 40, Tsinghua University
#
# The Haskell Implementation is as follows:
#
# ah' :: (RealFloat a) => Int -> Int -> [a] -> a -> a
# ah' i n xs x
#     | i == 0 = if x < xs !! 0 then 0 else (if x <= xs !! 1 then aleft (xs !! 0) (xs !! 1) x else 0)
#     | i == n = if x < xs !! (n - 1) then 0 else (if x <= xs !! n then aright (xs !! (n-1)) (xs !! n) x else 0)
#     | otherwise = if x < (xs !! (i - 1)) then 0 else (if x <= (xs !! i) then aright (xs !! (i - 1)) (xs !! i) x else (if x <= (xs !! i + 1) then aleft (xs !! i) (xs !! (i + 1)) x else 0))
#
# i = 0, 1, ... n (for a total of n+1 nodes)
# Indexes are mathematical since the mathematical definition starts with zero.
def aI(x, i, nodes):
    n = len(nodes) - 1
    if(n < 1):
        raise ValueError("There must be at least 2 nodes for Hermite Extrapolation (or any, actually.)")

    if(i == 0):
        return 0 if x < nodes[0] else (aL(a=nodes[0],b=nodes[1],x=x) if x < nodes[1] else 0)
    elif(i == n):
        return 0 if x < nodes[n-1] else (aR(a=nodes[n-1],b=nodes[n],x=x) if x <= nodes[n] else 0)
    else:
        return 0 if x < nodes[i-1] else (aR(a=nodes[i-1],b=nodes[i],x=x) if x < nodes[i] else (aL(a=nodes[i],b=nodes[i+1],x=x) if x < nodes[i+1] else 0)) 

def bI(x, i, nodes):
    n = len(nodes) - 1
    if(n < 1):
        raise ValueError("There must be at least 2 nodes for Hermite Extrapolation (or any, actually.)")

    if(i == 0):
        return 0 if x < nodes[0] else (bL(a=nodes[0],b=nodes[1],x=x) if x < nodes[1] else 0)
    elif(i == n):
        return 0 if x < nodes[n-1] else (bR(a=nodes[n-1],b=nodes[n],x=x) if x <= nodes[n] else 0)
    else:
        return 0 if x < nodes[i-1] else (bR(a=nodes[i-1],b=nodes[i],x=x) if x < nodes[i] else (bL(a=nodes[i],b=nodes[i+1],x=x) if x < nodes[i+1] else 0))

# Generate Equations
# Given xs, ys, generate the matrices to be solved in order to get the appropriate m-factors
# in the Splined Hermite representation Sh(x).
#
# Note: This is *not* very efficient as it stores a whole matrix.
# But it is generalized enough so we can abstract a solver out of the equation generator,
# so its for the greater good.
#
# Arrays xs and ys are mathematical indexes i=0..n
def genSplineMatrix(xs, ys):
    n = len(xs) - 1
    if(n != len(ys) - 1):
        raise ValueError("You must provide as many y-points as x-points to genSplineMatrix")

    M = mp.diag_(n+1, [2] * (n+1))

    # arrays lambdas & mus are mathematical indexes i=0..n
    lambdas = [None] * (n+1)
    lambdas[0] = 1
    lambdas[n] = 0
    M.setElem_(1, 2, lambdas[0])
    M.setElem_(n+1, n, 1)
    # h_i = x_i+1 - x_i

    mus = [None] * (n+1)
    mus[0] = 3 * (ys[1] - ys[0]) / (xs[1] - xs[0])
    mus[n] = 3 * (ys[n] - ys[n-1]) / (xs[n] - xs[n-1])

    for i in range(1, n):
        lambdas[i] = (xs[i] - xs[i-1]) / (xs[i+1] - xs[i-1])
        mus[i]     = 3 * ((ys[i] - ys[i-1]) * (1 - lambdas[i]) / (xs[i] - xs[i-1]) + (ys[i+1] - ys[i]) * lambdas[i] / (xs[i+1] - xs[i]))
        # and commit using the side-effect setElem_ function
        M.setElem_(i+1,i+2,lambdas[i])
        M.setElem_(i+1,i,1-lambdas[i])

    B = mp.fromList_(n+1,1,mus)

    return (M, B)

# Reticulate Splines (Actual Code)
# Generate Equations, Solve Equations, Construct Functions, and returns a Lambda
# representing the splined function.
#
# Luckily, Python does not do lazy evaluation, meaning that this lambda will be speedy
# instead of being tri-solved every call. :)
def spline_(xs, ys, x):
    gendMatrix = genSplineMatrix(xs, ys)
    mSolution  = solveTri(gendMatrix[0], gendMatrix[1]).toList() # ms...

    # To do a sum, we need parrying & foldl which Python does not do very well.
    # Hence we use functools.partial to do this in a tricky way,
    # making this spline_ function do the dirty work and wrapping it in a cleaner
    # spline function.
    # def aI(x, i, nodes) / bI
    res = 0
    for i, y in enumerate(ys):
        res += y * aI(x, i, xs) + mSolution[i] * bI(x, i, xs)

    return res

def spline(xs, ys):
    spp = partial(spline_, xs=xs, ys=ys)
    return lambda x: spp(x=x)

ex8 = spline([0,3,5,7,9,11,12,13,14,15], [0,1.2,1.7,2.0,2.1,2.0,1.8,1.2,1.0,1.6])

# For plotting the result
import matplotlib.pyplot as plt 
fig, ax = plt.subplots()
xs = [x/100 for x in range(1,1500)]
ys = [ex8(x) for x in xs]
ax.plot(xs, ys)
plt.show()

# For 0.1x variances
xs = [0,3,5,7,9,11,12,13,14,15]
for i in range(len(xs)):
    print("x =", xs[i] - 0.1, " y =", ex8(xs[i] - 0.1))
    print("x =", xs[i], " y =", ex8(xs[i]))
    print("x =", xs[i] + 0.1, " y =", ex8(xs[i] + 0.1))