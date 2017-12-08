####################################################################
# Computational Physics, 2017-18 Sem1
# HW-2 Ex-3
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
# Now Playing: The Scientist - Coldplay
#
# 2017.12.05
# Happy Birthday, Joy
####################################################################

debugLevel = 200

# Provide a seed
LCGSeed = None

# If no seed, then seed with time
if(LCGSeed == None):
    import time
    LCGSeed = int(time.time())


# LCGRand()
# Generate random number in [0, 1] using the "Linear Congruence Generator"
# (1951, D.H. Lehmer)
#
# Using Parameters from Numerical Recipes:
# m = 2**32, a = 1664525, c = 1013904223
# But this is customizable.
#
# It is IMPERATIVE to seed LCGRand() with a truly random seed value,
# before starting its usage.
# For the sole purposes of HW-3, LCG is seeded with the system time,
# but you SHOULD NOT USE THIS FOR CRYPTOGRAPHIC PURPOSES
# (hell, you should not be using RNGs at all!)

# Internal State.
LCGInternalState = LCGSeed
def LCGRand():
    global LCGSeed, LCGInternalState

    # Numerical Recipes Parameters for LCG
    a = 1664525
    c = 1013904223
    m = 2**32

    LCGInternalState = (a * LCGInternalState + c) % m
    return LCGInternalState / m

# Do not remove.
# This runs the LCGRand a few times to remove the initial bias
# caused by time being a not-so-fast-moving variable.
for i in range(1, 255):
    LCGRand()

pi = 3.141592653589793
e  = 2.718281828459045
# Utility: sine function (using complex repr. for speed)
def sin(theta):
    return 0.5*(e**(1j*theta) - e**(-1j*theta))/1j

def cos(theta):
    return 0.5*(e**(1j*theta) + e**(-1j*theta))

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
def solveSecant(f, a, b, xmin = None, xmax = None, maxIterations = 2000, maxE = 1e-15, threshold = 1e-4, isComplexAware = False):
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
        if(debugLevel >= 100):
            print("* diag solveSecant: stopped after", iterCount, " iterations to reach", currX, " res:", f(currX), flush=True)
        return None

    if(debugLevel >= 1000):
        print("* diag solveSecant: called", iterCount, " iterations to reach", currX, flush=True)
    
    return currX

# Ex-3 Specific Code below
# Ex3-(2)
# 1/4 each direction (l/r/u/d)
X = 0
Y = 0
print("** Ex3-(2)")
for i in range(50):
    rand = LCGRand()
    if(rand <= 0.25):
        X += 1
    elif(rand <= 0.5):
        X += -1
    elif(rand <= 0.75):
        Y += 1
    else:
        Y += -1
    print("(X, Y) = (" + X.__str__() + ", " + Y.__str__() + ") R = " + ((X**2 + Y**2)**(1/2)).__str__())

print("\n")

# Ex3-(2)b
# Randomized for each direction. Mapping [0, 1] -> [0, 2*pi] for the direction,
# then just multiply by sine and cosine

X = 0
Y = 0
print("** Ex3-(2)b")
for i in range(50):
    randT = LCGRand() * 2 * pi
    X += cos(randT).real
    Y += sin(randT).real
    print("(X, Y) = (" + X.__str__() + ", " + Y.__str__() + ") R = " + ((X**2 + Y**2)**(1/2)).__str__())

print("\n")

# Ex3-(3)
# 100 steps, output R and expectation value of R (E(R))
print("** Ex3-(2)a")
for i in range(100):
    randT = LCGRand() * 2 * pi
    X += cos(randT).real
    Y += sin(randT).real
    print("Step", i + 1, "R = " + ((X**2 + Y**2)**(1/2)).__str__())

print("\n")

# Perform NumberOfExperiments amount of experiments, accumulating-R to calculate Exp. Value
NumberOfExperiments = 30000
Rs = [0] * NumberOfExperiments

print("** Ex3-(3) NumberOfExperiments =", NumberOfExperiments, " wait...")
for i in range(NumberOfExperiments):
	X = 0
	Y = 0
	for i in range(100):
	    randT = LCGRand() * 2 * pi
	    X += cos(randT).real
	    Y += sin(randT).real
	    R = (X**2 + Y**2)**(1/2)
	    Rs[i] += R
print("** Result: ")
for i in range(100):
	Rs[i] = Rs[i] / NumberOfExperiments
	print(i+1, Rs[i])

print("\n")
print("* fitting ...")
# precalculate sum of -y_i * sqrt(n_i)
b = sum([Rs[i+1] * (i+1)**(1/2) for i in range(len(Rs) - 1)])
f = lambda a: sum(a * (i + 1) for i in range(len(Rs) - 1)) - b
ar = solveSecant(f, 0, 0.1)
print(ar, f(ar))

# Ex3-(3)b
# Internally calculate E(R) for a wide range of values and perform fitting on data.
# Rsum = 0
# X = 0
# Y = 0
# Rexps = [0] # padd so we have real mathematical indexes
# print("** Ex3-(3)b generating data...")
# for i in range(1, 500000):
#     randT = LCGRand() * 2 * pi
#     X += cos(randT).real
#     Y += sin(randT).real
#     R = (X**2 + Y**2)**(1/2)
#     Rsum += R
#     Rexp = (Rsum/(i))
#     Rexps.append(Rexp)
#     # print(i, Rexp)
# print("* fitting ...")
# # precalculate sum of -y_i * sqrt(n_i)
# b = sum([Rexps[i+1] * (i+1)**(1/2) for i in range(len(Rexps) - 1)])
# f = lambda a: sum(a * (i + 1) for i in range(len(Rexps) - 1)) - b
# ar = solveSecant(f, 0, 0.1)
# print(ar, f(ar))