
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

# DIAGNOSTICS:
#   Set to 0 for quiet (production mode)
#   If > 0, various levels (100, 200, 1000) will return gradually
#   verbose diagnostic messages for the purpose of debugging.
# Not all functions support a diagnostic output.
debugLevel = 0

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
e  = 2.718281828459045
# Utility: Cosine Function
# Uses Taylor Series Approximation for cosine
# Depends on:
#    factorial
# FIXME: 64 items is also an arbitrary precision limit.
# Improved 10-23-2017: If there is a modulo 2pi, remove it for higher
#   precision.
def cos(theta):
    if(isinstance(theta, complex)): # is complex, use complex definition
        return 0.5*(e**(1j*theta) + e**(-1j*theta))

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
    if(isinstance(theta, complex)): # is complex, use complex definition
        return 0.5*(e**(1j*theta) - e**(-1j*theta))/1j

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

# Trapezoid Formula for Integration
# Num integTrapezoid(f lambda x: y, a, b, h)
def integTrapezoid(f, a, b, n):
    global debugLevel
    
    result = 0
    h  = (b - a)/n
    
    # Sum over k = 1, ... n-1 for f(x_k)
    result = h/2 * (f(a) + f(b) + 2*sum([f(a+i*h) for i in range(1,n)]))

    if(debugLevel >= 200):
        print("* integTrapezoid: integrated from a=", a, "b=", b, " with n=", n, "partitions I=", result, flush=True)

    return result

####################################################################
# Subpoint (1)
# Solve saddle point equation for given p (pz) ts1 ~ ts6

# double|complex solveSecant(
#   lambda f: x, a, b, xmin = None, xmax = None, 
#   maxIterations = 1000, maxE = 1e-11,
#   threshold = 1e-5,
#   isComplexAware = False
# )
# Solves for roots with secant method with given initial values a, b
# With maximum iteration number maxIterations, or when deltaE <= maxE,
# or when x overflows beyond [xmin, xmax] closed range.
# Only returns a root if it does |f(x)| <= threshold, otherwise None.
# This prevents false positives.
#
# If isComplexAware is set to True, the search will also be Complex number aware.
# In this case, xmin and xmax restrictions are applied on real & imag parts
# respectively.
def solveSecant(f, a, b, xmin = None, xmax = None, maxIterations = 2000, maxE = 1e-15, threshold = 1e-7, isComplexAware = False):
    global debugLevel
    iterCount = 1
    currX = b
    lastX = a
    while(iterCount < maxIterations):
        if(xmin != None and xmax != None):
            if((isComplexAware and (currX.real < xmin.real and currX.imag < xmin.imag) or (currX.real > xmax.real and currX.imag > xmax.real)) or (not isComplexAware and (currX < xmin or currX > xmax))):
                if(abs(f(lastX)) > threshold):
                    return None
                else:
                    return lastX

        try:
            if(f(currX) - f(lastX) == 0): # underflow
                break

            (lastX, currX) = (currX, currX - (currX - lastX)/(f(currX) - f(lastX)) * f(currX))
        except OverflowError:
            if(debugLevel > 0):
                print("* diag solveSecant: **OVERFLOW ERROR**", flush=True)
            return None

        if(abs(lastX - currX) <= maxE):
            break
        iterCount += 1

    if(abs(f(currX)) > threshold or currX != currX): # last check for isNaN
        return None

    if(debugLevel >= 1000):
        print("* diag solveSecant: called", iterCount, " iterations to reach", currX, flush=True)
    
    return currX

# [double|complex] solveAllSecant(
#    lambda f: x, a, b, seedDistance = 0.01, maxE = 1e-11,
#    nudgeFactor = None,
#    isMP = True
# )
# Attempts to solve for ALL existing roots for (a mostly well-behaved) function
# By doing some very un-scientific seeding across [a, b]
# The seeds are spaced by seedDistance and initial values are spaced by a nudgeFactor
# within seedDistance, by default 1e4 of maxE
# Returns a list of all roots that are found, considering all roots within the distance
# of maxE the same root.
#
# If isComplexAware is specified, then Complex-space searching is used instead.
# The method starts from the reals then seeding across imX in [ia, ib] (specify or arbitrary
# defaults are used)
#
# If isMP is specified, multiprocessing.dummy-based Multi-Core processing will be used.
# Please specify the global variable mp_max_threads for the number of pools.
def solveAllSecant(f, a, b, seedDistance = 0.01, maxE = 1e-15, nudgeFactor = 1e4, isComplexAware = False, ia = -100, ib = 100):
    xs = [a + i * seedDistance for i in range(int((b - a)//seedDistance + 1))]
    if(isComplexAware):
        xs = [x + ia * 1j + i * seedDistance * 1j for i in range(int((ib-ia)//seedDistance + 1)) for x in xs]

    rs = [solveSecant(f, x, x + maxE * nudgeFactor, xmin = (a if not isComplexAware else a + ia * 1j), xmax = (b if not isComplexAware else b + ib * 1j), maxE = maxE, isComplexAware = isComplexAware) for x in xs]
    
    # comb through the results and retrieves maxE
    # also filtering through some results which don't make much sense
    rs = [r for r in rs if not r == None]

    # perform a filter which removes the errors by "rounding" the parts
    # don't learn from me... this is terrible I won't admit I wrote this
    targetPrec = abs(int(maxE.__str__().split("e")[1]))
    rs = [round(r.real, targetPrec) + round(r.imag, targetPrec) * 1j for r in rs]

    # eliminate trivially "same" solutions
    rs = list(set(rs))

    # if(not isComplexAware):
    #     rs.sort()
    # else:
    #     rs.sort(key = lambda z: z.real)
    #     rs.sort(key = lambda z: z.imag)
    #     # x < y then -1, x == y then 0, x > y then 1 is the sort definition
    #     # warning: this complex "sort" is useful only for eliminating duplicates.
    #     # it DOES NOT form a proper mathematical sort for complex numbers.
    
    # If you want to do maxE-circle based validation then change to len(rs)-1
    for i in range(len(rs)):
         if(isComplexAware):
             if(rs[i].real < a or rs[i].real > b or rs[i].imag < ia or rs[i].imag > ib):
                 rs[i] = None

    #     if(rs[i] != None and abs(rs[i] - rs[i+1]) <= maxE):
    #         # del rs[i]
    #         rs[i] = None

    # pass
    rs = [r for r in rs if not r == None]

    # diagnostics: remove in production
    if(debugLevel >= 200):
        print("*** solveAllSecant diag: we found", len(rs), "roots")
        for i in range(len(rs)):
            print("* x =", rs[i], " |f(x)|=", abs(f(rs[i])))

    return rs

# Ex-1 Constants in a.u.
A0 = 2.6509412245864
omega = 0.02847709533225
t0 = 0
tT = 441.27992928154773
Ip = 0.5

########################################################
# Actual code for Ex1-(1)
# Use x = omega * t instead as a variable, with range 0 ~ 12.566370614359173
Td1S = lambda t: 0.5 * (1 + A0 * sin(t * omega) * sin(t * omega / 4) * sin(t * omega / 4))**2 + Ip

# Omit the solutions where the Im part < 0
if __name__ == '__main__':
    print("* Exercise 1-(1):")
    print(solveAllSecant(f = Td1S, a = t0, b = tT, seedDistance = 50, maxE = 1e-12, nudgeFactor = 1e3, isComplexAware = True, ia = -0.1))
    print("\n\n")


########################################################
# Actual code for Ex1-(2)
# Abstracted for every Pz, where Ek=Pz^2/2
# Returns |Mp0|^2
def ex12(Pz):
    global A0, omega, t0, tT, Ip
    Ek = Pz**2 / 2
    Td1S = lambda t: 0.5 * (Pz + A0 * sin(t * omega) * sin(t * omega / 4) * sin(t * omega / 4))**2 + 0.5
    ts = solveAllSecant(f = Td1S, a = t0, b = tT, seedDistance = 50, maxE = 1e-9, nudgeFactor = 1e3, isComplexAware = True, ia = -0.1)

    TdTIS = lambda x: 0.5*(Pz + A0*sin(omega*x)*sin(omega*x/4)*sin(omega*x/4))**2 + Ip
    TS = lambda t: integTrapezoid(TdTIS, 0, t, 1000)

    Tid2ES = lambda t: 1/(A0*A0*omega*sin(omega*t)*(sin(omega*t/4))**4 + 0.5*A0*A0*omega*(sin(omega*t))**2*(sin(omega*t/4))**3*cos(omega*t/4)) * e**(1j * TS(t))

    fts = [Tid2ES(t) for t in ts]
    Mp0_spm = -(2 * Ip)**(5/4) / 2**(1/2) * sum(fts)

    return abs(Mp0_spm)**2

if __name__ == '__main__':
    print("* Exercise 1-(2):")
    from multiprocessing import pool

    partitions = 1000
    Pzs  = [0 + (2 / partitions) * i for i in range(partitions + 1)]
    Eks  = [Pz**2 / 2 for Pz in Pzs]

    MPool = pool.Pool(4)
    Mp0sqs = MPool.map(ex12, Eks)
    # Mp0sqs = [ex12(Pz) for Pz in Pzs]

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.plot(Eks, Mp0sqs)
    plt.show()
