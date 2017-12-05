####################################################################
# Computational Physics, 2017-18 Sem1
# HW-2 Ex-2
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
# This program is Python 3 (3.6.2 on MSYS_NT)
# compatible, and does not support Python 2.
#
# Now Playing: Kids In The Dark - All Time Low
####################################################################

from ArbitraryPrecision.Core import ArbitraryPrecision

# DIAGNOSTICS:
#   Set to 0 for quiet (production mode)
#   If > 0, various levels (100, 200, 1000) will return gradually
#   verbose diagnostic messages for the purpose of debugging.
# Not all functions support a diagnostic output.
debugLevel = 0

# HermitePoly(Int n, Numeric x)
# Computes Hermite Polynomial for given n, x. Is ArbitraryPrecision aware.
def HermitePoly(n, x):
    global debugLevel
    if(n == 0):
        return 1
    if(n == 1):
        return x*2
    if(n == 2):
        return x*x*4 - 2

    Hs = [1, x*2, x*x*4 - 2]
    for i in range(3, n + 1):
        Hs.append(x*2*Hs[i-1] - Hs[i-2]*2*(i - 1))

    if(debugLevel >= 1000):
        print("* diag HermitePoly: called for n=", n, "x=", x, flush=True)

    return Hs[n]

# dHermitePoly(Int n, Numeric x)
# Computes Derivative of nth Hermite Polynomial for given n, x.
# Is ArbitraryPrecision aware.
def dHermitePoly(n, x):
    return HermitePoly(n - 1, x) * (2 * n)

# double|ArbitraryPrecision newtonSolve(
#   lambda f: x, 
#   lambda df: x,
#   x0, xmin = None, xmax = None, 
#   maxIterations = 1000,
#   threshold = 1e-5
# )
# Solves for roots with Newton's method with given initial value x0
# With maximum iteration number maxIterations,
# or when x overflows beyond [xmin, xmax] closed range.
# Only returns a root if it does |f(x)| <= threshold, otherwise None.
# This prevents false positives.
#
# This function is ArbitraryPrecision aware,
# but does not support complex types in Python.
def newtonSolve(f, df, x0, xmin = None, xmax = None, maxIterations = 50, threshold = 1e-7):
    global debugLevel
    iterCount = 1
    currX = x0
    lastX = x0
    while(iterCount < maxIterations):
        if(xmin != None and xmax != None):
            if(currX < xmin or currX > xmax):
                if(abs(f(lastX)) > threshold):
                    return None
                else:
                    return lastX

        dfx = df(currX)
        if(dfx == 0):
            if(abs(f(lastX)) > threshold):
                return None
            else:
                return lastX

        (lastX, currX) = (currX, currX - f(currX)/df(currX))
        iterCount += 1

        print("*", currX, flush=True)

    if(abs(f(currX)) > threshold or currX != currX): # last check for isNaN
        return None

    if(debugLevel >= 200):
        print("* diag newtonSolve: called", iterCount, " iterations to reach", currX, flush=True)
    
    return currX

print(newtonSolve(lambda x: HermitePoly(288, x), lambda x: dHermitePoly(288, x), ArbitraryPrecision(12)))

# print(HermitePoly(1500, ArbitraryPrecision(5.122)))
