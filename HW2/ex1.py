####################################################################
# Computational Physics, 2017-18 Sem1
# HW-2 Ex-1
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
# Now Playing: Learn to Love Again - Lawson
####################################################################

####################################################################
# Utility Functions
####################################################################
# Utility: Factorial Function
# factorial(n) returns n! (for n in natural number)
# With Caching w/ Automatic Compilation (to avoid RecursionError)
factorial_cache = {0: 1, 1: 1}
def factorial(n, isCache = None):
    if(n == 0):
        return 1

    if(n > 500 and isCache == None):
        for i in range(249, n//2):
            factorial(i*2, True)

    cache = factorial_cache.get(n, None)
    if(cache == None):
        factorial_cache[n] = n * factorial(n - 1, isCache)
    return factorial_cache[n]

pi = 3.141592653589793
# Utility: Cosine Function
# Uses Taylor Series Approximation for cosine
# Depends on:
#    factorial
# FIXME: 64 items is also an arbitrary precision limit.
# Improved 10-23-2017: If there is a modulo 2pi, remove it for higher
#   precision.
def cos(theta):
    result = 1
    theta  = abs(theta) # symmetry helps
    while(theta > 2*pi):
        theta += (-1)*2*pi

    if(abs(theta-pi/2) < 1e-11 or abs(theta - 3*pi/2) < 1e-11):
        return 0
    elif(abs(theta-pi) < 1e-11):
        return (-1)
    elif(abs(theta) < 1e-11 or abs(theta-2*pi) < 1e-11):
        return 1

    if(theta > pi):
        return (-1) * cos(theta - pi) # reduce...
    if(theta > pi/2):
        return (-1) * sin(theta - pi/2)

    for i in range(1, 64):
        result += (-1)**i * theta**(i*2) / factorial(i*2)
    return result

# Utility: Sine Function
# Uses Taylor Series Approximation for cosine
# Depends on:
#    factorial
# FIXME: 64 items is also an arbitrary precision limit.
# Improved 10-23-2017: If there is a modulo 2pi, remove it for higher
#   precision.
def sin(theta):
    result = 0
    if(theta < 0):
        while(theta < (-1)*2*pi):
            theta += 2*pi 
    else:
        while(theta > 2*pi):
            theta += (-1)*2*pi

    for i in range(1, 64):
        result += (-1)**(i+1) * theta**(i*2-1) / factorial(i*2-1)
    return result

####################################################################
# Subpoint (1)
# Solve saddle point equation for given p (pz) ts1 ~ ts6

# double|complex solveSecant(
#   lambda f: x, a, b, xmin = None, xmax = None, 
#   maxIterations = 1000, maxE = 1e-11,
#   threshold = 1e-5
# )
# Solves for roots with secant method with given initial values a, b
# With maximum iteration number maxIterations, or when deltaE <= maxE,
# or when x overflows beyond [xmin, xmax] closed range.
# Only returns a root if it does |f(x)| <= threshold, otherwise None.
# This prevents false positives.
def solveSecant(f, a, b, xmin = None, xmax = None, maxIterations = 2000, maxE = 1e-15, threshold = 1e-5):
    iterCount = 1
    currX = b
    lastX = a
    while(iterCount < maxIterations):
        if(xmin != None and xmax != None and (currX < xmin or currX > xmax)):
            if(abs(f(lastX)) > threshold):
                return None
            else:
                return lastX
        if(f(currX) - f(lastX) == 0): # underflow
            break
        (lastX, currX) = (currX, currX - (currX - lastX)/(f(currX) - f(lastX)) * f(currX))
        if(abs(lastX - currX) <= maxE):
            break
        iterCount += 1

    if(abs(f(currX)) > threshold):
        return None

    # print("* diag solveSecant: called", iterCount, " iterations", flush=True)
    return currX

# [double|complex] solveAllSecant(
#    lambda f: x, a, b, seedDistance = 0.01, maxE = 1e-11,
#    nudgeFactor = None
# )
# Attempts to solve for ALL existing roots for (a mostly well-behaved) function
# By doing some very un-scientific seeding across [a, b]
# The seeds are spaced by seedDistance and initial values are spaced by a nudgeFactor
# within seedDistance, by default 1e4 of maxE
# Returns a list of all roots that are found, considering all roots within the distance
# of maxE the same root.
def solveAllSecant(f, a, b, seedDistance = 0.01, maxE = 1e-15, nudgeFactor = 1e4):
    xs = [a + i * seedDistance for i in range(int((b - a)//seedDistance + 1))]
    rs = [solveSecant(f, x, x + maxE * nudgeFactor, xmin = a, xmax = b, maxE = maxE) for x in xs]
    
    # comb through the results and retrieves maxE
    # also filtering through some results which don't make much sense
    rs = [r for r in rs if not r == None]
    rs.sort()
    for i in range(len(rs) - 1):
        if(rs[i] != None and abs(rs[i] - rs[i+1]) <= maxE):
            # del rs[i]
            rs[i] = None

    return [r for r in rs if not r == None]

# Ex-1 Constants.
A0 = 2.6509412245864
omega = 0.02847709533225
t0 = 0
tT = 441.27992928154773

# Actual code for Ex1-(1)
# d1S = 0.5 * Txd1S**2 = 0 equiv. Txd1S == 0
# Use x = omega * t instead as a variable, with range 0 ~ 12.566370614359173
Txd1S = lambda x: 1 + A0 * sin(x) * sin(x / 4) * sin(x / 4)
print(solveAllSecant(Txd1S, omega * t0, omega * tT, 1e-5, maxE = 1e-10, nudgeFactor = 1e3))