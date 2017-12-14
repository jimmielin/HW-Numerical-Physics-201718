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
debugLevel = 201

# HermitePoly(Int n, Numeric x)
# Computes Hermite Polynomial for given n, x. Is ArbitraryPrecision aware.
#
# If WithDeriv set to True, returns a tuple instead with (H[x], H'[x]) at degree-n.
# (useful for Newton Solving)
def HermitePoly(n, x, WithDeriv = False):
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

    if(WithDeriv):
        return (Hs[n], Hs[n-1]*(2*n))
    else:
        return Hs[n]

# dHermitePoly(Int n, Numeric x)
# Computes Derivative of nth Hermite Polynomial for given n, x.
# Is ArbitraryPrecision aware.
def dHermitePoly(n, x):
    return HermitePoly(n - 1, x) * (2 * n)

# double|ArbitraryPrecision bisectionSolve(
#   lambda f: x,
#   xmin = None, xmax = None,
#   maxIterations = 1000,
#   maxE = 1e-10,
#   threshold = 1e-5
#)
#
# A bisection-solver with a behavior similarly
# compatible with newtonSolve, with exception of x0 omitted
def bisectionSolve(f, xmin = None, xmax = None, maxIterations = 500, maxE = 1e-10, threshold = 1e-7):
    global debugLevel
    iterCount = 1
    a = xmin
    b = xmax
    if(f(a).sgn() == f(b).sgn()):
        return None

    while(iterCount < maxIterations):
        c = (a + b)/2
        fa = f(a)
        fc = f(c)
        if(fa.sgn() == fc.sgn()):
            a = c
        else:
            b = c
        iterCount += 1

    fxx = f(a)
    if(threshold != None and abs(fxx) > threshold or a != a): # last check for isNaN
        if(debugLevel >= 100):
            print("* diag bisectionSolve: after", iterCount, " eliminated", a, " f(x)=", fxx, flush=True)
        return None

    if(debugLevel >= 200):
        print("* diag bisectionSolve: called", iterCount, " iterations to reach", a, flush=True)
    
    return a



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
#
# If df = None, f is assumed to be a function that returns a tuple with (f, f') instead.
# This may be useful for optimizations.
#
# If threshold = None, which is useful for some VERY EXTREME oscillating functions,
# then threshold calculation is completely ignored & all results are returned. Beware!
def newtonSolve(f, df, x0, xmin = None, xmax = None, maxIterations = 50, threshold = 1e-7):
    global debugLevel
    iterCount = 1
    currX = x0
    lastX = x0
    while(iterCount < maxIterations):
        flx = (abs(f(lastX)) if df != None else abs(f(lastX)[0]))
        if(xmin != None and xmax != None):
            if(currX < xmin or currX > xmax):
                if(threshold != None and flx > threshold):
                    return None
                else:
                    return lastX

        if(df != None):
            dfx = df(currX)
            fx  = f(currX)
        else:
            fc  = f(currX)
            fx  = fc[0]
            dfx = fc[1]

        if(dfx == 0):
            if(threshold != None and flx > threshold):
                return None
            else:
                return lastX

        (lastX, currX) = (currX, currX - fx/dfx)
        iterCount += 1

    fxx = (f(currX) if df != None else f(currX)[0])
    if(threshold != None and abs(fxx) > threshold or currX != currX): # last check for isNaN
        if(debugLevel >= 100):
            print("* diag newtonSolvd: after", iterCount, " eliminated", currX, " f(x)=", fxx, flush=True)
        return None

    if(debugLevel >= 200):
        print("* diag newtonSolve: called", iterCount, " iterations to reach", currX, flush=True)
    
    return currX

# Multithreading Caller
# This is just a wrapper to avoid functools.partial, solely designed for solveAllNewton on the Hermite Poly
def _hermite_solve_mpCaller_288(x):
    return newtonSolve(f = lambda x: HermitePoly(288, x, True), df = None, x0 = x, xmin = ArbitraryPrecision(1, -15, InternalAware = True), xmax = ArbitraryPrecision(23.378), maxIterations = 20, threshold = None)


# [double|ArbitraryPrecision] solveAllNewton(
#    lambda f: x,
#    lambda df: x,
#    a, b, partitions = 1000, maxE = 1e-15
# )
# Attempts to solve for ALL existing roots for (a mostly well-behaved) function
# By doing even-spaced seeding across [a, b] in Ordered Numeric set.
#
# The seeds are spaced evenly into partitions.
# Returns a list of all roots that are found, considering all roots within the distance
# of maxE the same root. Only the exponential part is considered.
#
# This function is ArbitraryPrecision aware.
def solveAllNewton(f, df, a, b, partitions = 1000, maxE = 1e-5):
    seedDistance = (b - a)/partitions
    xs = [a + seedDistance * i for i in range(partitions)]

    # rs = [newtonSolve(f, df, x0 = x, xmin = a, xmax = b, maxIterations = 5, threshold = None) for x in xs]
    ## THIS IS FOR HERMITE SOLVING ONLY ##
    from multiprocessing import pool
    MPool = pool.Pool(8)
    rs = MPool.map(_hermite_solve_mpCaller_288, xs)
    
    # comb through the results and retrieves maxE
    # also filtering through some results which don't make much sense
    rs = [r for r in rs if not r == None]

    # eliminate trivially "same" solutions
    rs = list(set(rs))

    rs.sort()

    for i in range(len(rs) - 1):
        if(rs[i] != None and abs(rs[i] - rs[i+1]) <= maxE):
            # del rs[i]
            rs[i] = None

    # pass
    rs = [r for r in rs if not r == None]

    # diagnostics: remove in production
    if(debugLevel >= 200):
        print("*** solveAllNewton diag: we found", len(rs), "roots")
        for i in range(len(rs)):
            print("* x =", rs[i])

    return rs

# Multithreading Caller
# This is just a wrapper to avoid functools.partial, solely designed for solveAllBisection on the Hermite Poly with a fixed stepping size 0.1
def _hermite_solve_mpCaller_288_bis_1em1(x):
    return bisectionSolve(f = lambda x: HermitePoly(288, x), xmin = x, xmax = x + 0.1, maxIterations = 100, maxE = 1e-15, threshold = None)

def _hermite_solve_mpCaller_96_bis_1em1(x):
    return bisectionSolve(f = lambda x: HermitePoly(96, x), xmin = x, xmax = x + 0.1, maxIterations = 100, maxE = 1e-15, threshold = None)

def _hermite_solve_mpCaller_48_bis_1em1(x):
    return bisectionSolve(f = lambda x: HermitePoly(48, x), xmin = x, xmax = x + 0.1, maxIterations = 100, maxE = 1e-15, threshold = None)

def _hermite_solve_mpCaller_24_bis_1em1(x):
    return bisectionSolve(f = lambda x: HermitePoly(24, x), xmin = x, xmax = x + 0.1, maxIterations = 100, maxE = 1e-15, threshold = None)


# [double|ArbitraryPrecision] solveAllBisection(
#    lambda f: x,
#    lambda df: x,
#    a, b, distance = 0.01, partitions = (b-a)//distance + 1, maxE = 1e-15
# )
# Attempts to solve for ALL existing roots for (a mostly well-behaved) function
# By doing even-spaced seeding across [a, b] in Ordered Numeric set.
#
# The seeds are spaced evenly into partitions of distance distance.
# Returns a list of all roots that are found, considering all roots within the distance
# of maxE the same root. Only the exponential part is considered.
#
# This function is ArbitraryPrecision aware, BUT
# if you use ArbitraryPrecision you MUST provide BOTH distance and partitions
def solveAllBisection(f, df, a, b, distance, partitions, maxE = 1e-5):
    if(partitions == None):
        partitions = (b - a) // distance + 1
    xs = [a + distance * i for i in range(partitions)]

    # rs = [newtonSolve(f, df, x0 = x, xmin = a, xmax = b, maxIterations = 5, threshold = None) for x in xs]
    ## THIS IS FOR HERMITE SOLVING ONLY ##
    from multiprocessing import pool
    MPool = pool.Pool(8)
    rs = MPool.map(_hermite_solve_mpCaller_288_bis_1em1, xs)
    
    # comb through the results and retrieves maxE
    # also filtering through some results which don't make much sense
    rs = [r for r in rs if not r == None]

    # eliminate trivially "same" solutions
    rs = list(set(rs))

    rs.sort()

    for i in range(len(rs) - 1):
        if(rs[i] != None and abs(rs[i] - rs[i+1]) <= maxE):
            # del rs[i]
            rs[i] = None

    # pass
    rs = [r for r in rs if not r == None]

    # diagnostics: remove in production
    if(debugLevel >= 200):
        print("*** solveAllBisection diag: we found", len(rs), "roots")
        for i in range(len(rs)):
            print("* x =", rs[i])

    return rs

if __name__ == '__main__':
    #print(solveAllNewton(lambda x: HermitePoly(288, x, True), None, ArbitraryPrecision(1, -15, InternalAware = True), ArbitraryPrecision(23.4), 500))

    # print(solveAllBisection(lambda x: HermitePoly(288, x, True), None, ArbitraryPrecision(1, -6, InternalAware = True), ArbitraryPrecision(24), 0.1, 240))

    #print(bisectionSolve(lambda x: HermitePoly(288, x), ArbitraryPrecision(11), ArbitraryPrecision(30), threshold = None))

    print(HermitePoly(288, ArbitraryPrecision(5.122)))