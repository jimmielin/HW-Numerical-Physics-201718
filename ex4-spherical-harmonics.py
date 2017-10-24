####################################################################
# Computational Physics, 2017-18 Sem1
# HW-1 Ex-4
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
# Now Playing: Wait - M83
#              Beautiful Now - Zedd
####################################################################

####################################################################
# Numerical Differentiation
# Given a function f(x), get a numerical f'(x) result
# The timestep is dynamically assigned.
# f(x) must accept continuous values.
####################################################################

# Arbitrarily defined timestep
h_default = 0.0001

# 中心差商近似
# derivF_c(f = lambda x: f(x), x, h)
def derivF_c(f, x, h = None):
    if(h == None):
        h = h_default
    return (f(x + h) - f(x - h))/(2 * h)

# 二阶差商近似
# deriv2F_c(f = lambda x: f(x), x, h)
def deriv2F_c(f, x, h = None):
    if(h == None):
        h = h_default
    return (f(x + h) - 2 * f(x) + f(x - h))/(h * h)


# get f(level) derivative
# this is VERY inaccurate and should only be used for 1-3 levels
def derivF_c_cycle(f, x, level = 1, h = None):
    if(level == 0):
        return f(x)
    elif(level == 2):
        return deriv2F_c(f, x, h)
    elif(level == 1):
        return derivF_c(f, x, h)

    return derivF_c_cycle(lambda x: derivF_c(f, x), x, level - 1, h)

####################################################################
# Computing Legendre Polynomials
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

# Utility: Double Factorial Function
# doubleFactorial(n) returns n!! (double factorial) (for n in natural number)
# If acc (accumulator) is given, then the result is multiplied from the acc.
#  * WARNING: If you provide a acc, tail recursion is disabled (so is caching)
#      (doesn't really matter, Python's tail recursion is basically crap anyway)
doubleFactorial_cache = {0: 1, 1: 1, 2: 2}
def doubleFactorial(n, acc = 1, isCache = None):
    if(n == 0 or n == 1):
        return 1
    elif(n == 2):
        return 2

    if(acc != 1):
        while(n != 0 and n != 1):
            # print("* doubleFactorial diag: intermediate result = ", acc, flush=True)
            acc = acc * n
            n = n - 2
        return acc

    if(n > 250 and isCache == None):
        for i in range(124, n//2):
            doubleFactorial(i*2+n%2, isCache=True)

    cache = doubleFactorial_cache.get(n, None)
    if(cache == None):
        doubleFactorial_cache[n] = n * doubleFactorial(n - 2, isCache=isCache)
        #print("* doubleFactorial diag: cache missed n = ", n, flush=True)
    return doubleFactorial_cache[n]

# Utility: FactorialFromTo
# If from > to (this is an acceptable usage.), then the inverse is given and precisely calculated, float-first
# If acc (accumulator) is given, then the result is multiplied from the acc.
def factorialFromTo(a, b, acc = 1):
    result = acc
    if(a == b):
        return a

    if(a < b):
        for i in range(a + 1, b + 1):
            result = result * i
    elif(a > b):
        for i in range(b + 1, a + 1):
            result = result / i
            #print("* factorialFromTo diag: intermediate result at i = ", i, " result = ", result, flush=True)

    return result

# LegendrePolyNA (Non-Associated)
# With Caching w/ Pre-Compilation for large N (supports N up to 10000 as tested)
legendrePolyNA_cache = {0: {0: 1}, 1: {0: 0}}
def legendrePolyNA(n, x, isCache = None):
    global legendrePolyNA_cache
    # By definition (don't use this, it's completely not good)
    # return 1/(2**n * factorial(n)) * derivF_c_cycle(lambda x: (x**2 - 1)**n, x, n)

    # Using the recursive relation for Pn(x):
    # (n)P_n = (2n-1)xP_n-1 - (n-1)P_n-2, P0=1, P1=x
    # we can have better results..

    # If n is very large, maximum Recursion Depth will throw an error
    # so a cache should be built first, gradually
    if(n > 500 and isCache == None):
        for i in range(249, n//2):
            # print("* legendrePolyNA diag: building cache for n = ", n, isCache, flush=True)
            legendrePolyNA(i*2, x, True)

    if(n == 0):
        return 1
    elif(n == 1):
        return x

    cache = legendrePolyNA_cache.get(n, None)
    if(cache == None):
        legendrePolyNA_cache[n] = {}
        # print("* legendrePolyNA diag: cache missed for n = ", n, flush=True)
        legendrePolyNA_cache[n][x] = (2*n-1)/n*x*legendrePolyNA(n-1, x, isCache) - (n-1)/n*legendrePolyNA(n-2, x, isCache)
    else:
        cache = legendrePolyNA_cache[n].get(x, None)
        if(cache == None):
            # print("* legendrePolyNA diag: cache missed for n = ", n, flush=True)
            legendrePolyNA_cache[n][x] = (2*n-1)/n*x*legendrePolyNA(n-1, x, isCache) - (n-1)/n*legendrePolyNA(n-2, x, isCache)

    # print("* legendrePolyNA diag: n = ", n, ", x = ", x, flush=True)
    return legendrePolyNA_cache[n][x]

# LegendrePolyA (Associated Legendre Polynomial)
#    m
#  P    = ...
#    l
# With Caching w/ Automatic Pre-Compilation
# Caching structure: _cache[l][m][x]
#
# An improvement in 10-24-2017 uses Heuristic methods for choosing a recursion strategy, but it is not functional yet.
# Improved 10-24-2017: Supports m < 0 using factorialFromTo
legendrePolyA_cache = {}
def legendrePolyA(l, m, x, isCache = None):
    global legendrePolyA_cache

    # Define a few existing formulae we can calculate directly,
    # helpful for ending recursive results
    if(abs(m) > l):
        return 0

    if(m == 0):
        return legendrePolyNA(l, x)
    elif(m < 0):
        return (-1)**m * factorialFromTo(l-m, l+m, legendrePolyA(l, (-1)*m, x, isCache=isCache))

    if(l == m):
        return (-1)**l * doubleFactorial(2*l - 1) * ((1 - x**2)**(l / 2))
    #if(l-1 == m):
    #   return x*(2*l-1) * legendrePolyA(l-1, l-1, x, isCache)

    if(m == 1):
        if(l == 1):
            return (-1)*(1-x*x)**(1/2)
        if(l == 2):
            return (-1)*3*x*(1-x*x)**(1/2)

    # print("* legendrePolyA diag: l = ", l, ", m = ", m, ", x = ", x, ", isCache = ", isCache, flush=True)
    # Traditional Recursion:
    # This is a bit slow because it only decrements L
    # (2*l-1)/(l-m)*x*legendrePolyA(l-1, m, x) - (l+m-1)/(l-m)*legendrePolyA(l-2, m, x)
    # But it works for small L very nicely with higher precision,
    # also it is used for caching on large L
    # (2*l-1)/(l-m)*x*legendrePolyA(l-1, m, x) - (l+m-1)/(l-m)*legendrePolyA(l-2, m, x)

    # Pre-compiled Caching is used
    # to reduce recursion exponential depth growth.
    if(l > 25 and isCache == None):
        for i in range(24, l):
            # print("* legendrePolyA diag: building cache for l = ", i, " m = ", m, isCache, flush=True)
            legendrePolyA(i, m, x, True)


    # Another Recursion does incrementing & incrementing of M, while decrementing L
    # smartly, making use of legendryPolyNA,
    # used when M is sufficiently large.
    # and the fact that Plm=0 when |m|>l
    # ((1-x**2)**(1/2))*(-1)/(2*m)*(legendrePolyA(l-1, m+1, x)+(l+m-1)*(l+m)*legendrePolyA(l-1,m-1,x))

    # Heuristic for L/M-combo-values (if sufficiently large, use different strategy)
    #if(abs(l - m) < 3):
    #   # strategy = "dLiM" # decrement L, incrementing M *BRANCHES OUT EXPONENTIALLY
    #else:
    #   strategy = "dL_M" # decrement L only
    strategy = "dL_M"

    cache = legendrePolyA_cache.get(l, None)
    if(cache == None): # w/o l
        legendrePolyA_cache[l] = {}
        legendrePolyA_cache[l][m] = {}
        if(strategy == "dL_M"):
            legendrePolyA_cache[l][m][x] = (2*l-1)/(l-m)*x*legendrePolyA(l-1, m, x, isCache) - (l+m-1)/(l-m)*legendrePolyA(l-2, m, x, isCache)
        elif(strategy == "dLiM"):
            legendrePolyA_cache[l][m][x] = ((1-x**2)**(1/2))*(-1)/(2*m)*(legendrePolyA(l-1, m+1, x)+(l+m-1)*(l+m)*legendrePolyA(l-1,m-1,x))
        # print("* legendrePolyA diag: cache missed l = ", l, flush=True)
    else: # w l
        cache = legendrePolyA_cache[l].get(m, None)
        if(cache == None): # w/o m
            legendrePolyA_cache[l][m] = {}
            if(strategy == "dL_M"):
                legendrePolyA_cache[l][m][x] = (2*l-1)/(l-m)*x*legendrePolyA(l-1, m, x, isCache) - (l+m-1)/(l-m)*legendrePolyA(l-2, m, x, isCache)
            elif(strategy == "dLiM"):
                legendrePolyA_cache[l][m][x] = ((1-x**2)**(1/2))*(-1)/(2*m)*(legendrePolyA(l-1, m+1, x)+(l+m-1)*(l+m)*legendrePolyA(l-1,m-1,x))
            # print("* legendrePolyA diag: cache missed l = ", l, ", m = ", m, flush=True)
        else:
            cache = legendrePolyA_cache[l][m].get(x, None)
            if(cache == None): # w/o x
                if(strategy == "dL_M"):
                    legendrePolyA_cache[l][m][x] = (2*l-1)/(l-m)*x*legendrePolyA(l-1, m, x, isCache) - (l+m-1)/(l-m)*legendrePolyA(l-2, m, x, isCache)
                elif(strategy == "dLiM"):
                    legendrePolyA_cache[l][m][x] = ((1-x**2)**(1/2))*(-1)/(2*m)*(legendrePolyA(l-1, m+1, x)+(l+m-1)*(l+m)*legendrePolyA(l-1,m-1,x))
                # print("* legendrePolyA diag: cache missed l = ", l, ", m = ", m, ", x = ", x, ", isCache = ", isCache, flush=True)
            #else:
            #   print("* legendrePolyA diag: cache hit")

    # print("* legendrePolyNA diag: n = ", n, ", x = ", x, flush=True)
    return legendrePolyA_cache[l][m][x]

    # Given the identity:
    # (l-m)P(m,l) = (2l-1)xP(m,l-1) - (l+m-1)P(m,l-2)
    # We can lower the l-parameters until it is equal to m
    # return ((1-x**2)**(1/2))*(-1)/(2*m)*(legendrePolyA(l-1, m+1, x)+(l+m-1)*(l+m)*legendrePolyA(l-1,m-1,x))

    # By definition:
    # return (1 - x**2)**(m/2) * derivF_c_cycle(lambda x: legendrePolyNA(l, x), x, m)

####################################################################
# Computing Spherical Harmonics
####################################################################
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

# Utility: Exponential Function
# Uses Taylor Series Expansion
# Depends on:
#    factorial
# FIXME: 32 items an arbitrary precision limit
# Improved 10-23-2017: For purely imaginary inputs, the algorithm uses Euler.
def exp(z):
    if(z.imag != 0 and z.real == 0):
        return cos(z.imag) + 1j*sin(z.imag)
    else:
        result = 1
        for i in range(1, 32):
            result += z**i / factorial(i)
        return result

# _irSpFactorial(l)
# Special Heuristic Shim for M=L-1
# Supports large L, up to 1M
def _irSpFactorial(l):
    # (2l-3)!!(2l-3)!!/(2l-1)! = ((2l-3)(2l-5)...)^2/(2l-1)(2l-2)(2l-3)...
    # u*u/(u+2)
    result = 1
    ptr2 = 2*l - 3
    ptr1 = 2*l - 1
    while(ptr2 != -1):
        result = result * ptr2 / ptr1 
        if(ptr1 % 2 == 0):
            ptr2 = ptr2 - 2
        ptr1 = ptr1 - 1
    return result

# SphericalHarmonics(l, m, theta, phi)
# Computes spherical harmonics Ylm(theta, phi) for given parameters.
# Depends on the following functions:
#    legendrePolyA
#    cos (cosine function)
#    exp (exponential function)
#
# Notes:
#    ir prefix means "Intermediate Result" (Hungarian Notation)
def SphericalHarmonics(l, m, theta, phi):
    # Special Heuristic Shim for large M
    if(m == l-1):
        x = cos(theta)
        irRS = (2*l+1)/(4*pi)*_irSpFactorial(l)*x*x*(2*l-1)**2*(1-x**2)**(l-1)
        return (-1)**(l-1) * (-1 if x < 0 else 1) * (irRS)**(1/2) * exp(1j * m * phi)

    # factorial(l-m)/factorial(l+m) replace with 1/factorialFromTo(l-m, l+m) should be more accurate
    lPA = legendrePolyA(l, m, cos(theta))

    # power 1/2 of complex numbers / negative numbers is completely inaccurate in Python,
    # reaching absurd levels of error in the Real part (up to 1/2 of number)
    # Which is why we need a "safe", official, sanctioned square root.
    if(lPA < 0):
        lPA_ssqrt = ((-1)*lPA)**(1/2) * 1j
    else:
        lPA_ssqrt = lPA**(1/2)

    #print("* SphericalHarmonics diag: l=", l, " m=", m, " theta=", theta, " phi=", phi, flush=True)
    #print("* SphericalHarmonics diag: lPA intermediate result lPA=", lPA, " (sqrt)=", lPA_ssqrt, flush=True)

    #if(m > 50): # Heuristic for factorial explosion and custom accumulator
    #   ilr_factorialFromTo = factorialFromTo(l-m, l+m)
    #   if(ilr_factorialFromTo == 1e500): # +infty

    # Rely on some heuristics depending on intermediate results,
    # depending on what we obtain as the intermediate legendrePolyA result (e+150 is the maximum)
    if(abs(lPA) > 1e+150):
        irRS = factorialFromTo(l+m, l-m, (2*l+1) * lPA)/(4*pi)
        if(irRS < 0):
            irRS_ssqrt = ((-1)*irRS)**(1/2) * 1j
        else:
            irRS_ssqrt = irRS**(1/2)

        return irRS_ssqrt * lPA_ssqrt * exp(1j * m * phi)
    #elif(abs(lPA) < 1e-50 and m > 55): # Factorials grow large VERY fast and accumulators have to be used
    #   return (factorialFromTo(l+m, l-m, (2*l+1) * 1e300)/(4*pi))**(1/2) * lPA * 1e-300 * exp(1j * m * phi)
    else:
        irRS = factorialFromTo(l+m, l-m, (2*l+1) * lPA**2)/(4*pi)
        return (irRS)**(1/2) * (1 if lPA > 0 else (-1)) * exp(1j * m * phi)


#import time
#start_time = time.time()
#elapsed_time = time.time() - start_time
#print("Execution time is ", elapsed_time, " seconds.")

l = 1000

# print(SphericalHarmonics(l, 1, pi/1000, pi/6))
# print(SphericalHarmonics(l, 1, 3*pi/10, pi/6))
# print(SphericalHarmonics(l, 1, 501*pi/1000, pi/6))

# print(SphericalHarmonics(l, int(l/100), pi/1000, pi/6))
# print(SphericalHarmonics(l, int(l/100), 3*pi/10, pi/6))
# print(SphericalHarmonics(l, int(l/100), 501*pi/1000, pi/6))

# print(SphericalHarmonics(l, int(l/10), pi/1000, pi/6))
# print(SphericalHarmonics(l, int(l/10), 3*pi/10, pi/6))
# print(SphericalHarmonics(l, int(l/10), 501*pi/1000, pi/6))

print(SphericalHarmonics(l, l-1, pi/1000, pi/6))
print(SphericalHarmonics(l, l-1, 3*pi/10, pi/6))
print(SphericalHarmonics(l, l-1, 501*pi/1000, pi/6))
