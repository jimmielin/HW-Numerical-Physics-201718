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
def factorial(n):
	if(n == 0):
		return 1
	return n * factorial(n - 1)

# Utility: Double Factorial Function
# doubleFactorial(n) returns n!! (double factorial) (for n in natural number)
def doubleFactorial(n):
	if(n == 0 or n == 1):
		return 1
	elif(n == 2):
		return 2
	return n * doubleFactorial(n - 2)

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
		legendrePolyNA_cache[n][x] = (2*n-1)/n*x*legendrePolyNA(n-1, x) - (n-1)/n*legendrePolyNA(n-2, x)
	else:
		cache = legendrePolyNA_cache[n].get(x, None)
		if(cache == None):
			# print("* legendrePolyNA diag: cache missed for n = ", n, flush=True)
			legendrePolyNA_cache[n][x] = (2*n-1)/n*x*legendrePolyNA(n-1, x) - (n-1)/n*legendrePolyNA(n-2, x)

	# print("* legendrePolyNA diag: n = ", n, ", x = ", x, flush=True)
	return legendrePolyNA_cache[n][x]

# LegendrePolyA (Associated Legendre Polynomial)
#    m
#  P    = ...
#    l
# With Caching
# Caching structure: _cache[l][m][x]
legendrePolyA_cache = {}
def legendrePolyA(l, m, x):
	global legendrePolyA_cache

	# Define a few existing formulae we can calculate directly,
	# helpful for ending recursive results
	if(abs(m) > l):
		return 0

	if(m == 0):
		return legendrePolyNA(l, x)

	if(l == m):
		return (-1)**l * doubleFactorial(2*l - 1) * (1 - x**2)**(l / 2)
	if(l-1 == m):
		return x*(2*l-1) * legendrePolyA(l-1, l-1, x)

	if(m == 1):
		if(l == 1):
			return (-1)*(1-x*x)**(1/2)
		if(l == 2):
			return (-1)*3*x*(1-x*x)**(1/2)

	# print("* legendrePolyA diag: l = ", l, ", m = ", m, ", x = ", x, flush=True)
	# Old Recursion:
	# This is a bit slow because it only decrements L
	# (2*l-1)/(l-m)*x*legendrePolyA(l-1, m, x) - (l+m-1)/(l-m)*legendrePolyA(l-2, m, x)
	# But it works for small L very nicely with higher precision.

	# New Recursion does incrementing & incrementing of M, while decrementing L
	# smartly, making use of legendryPolyNA
	# and the fact that Plm=0 when |m|>l
	# ((1-x**2)**(1/2))*(-1)/(2*m)*(legendrePolyA(l-1, m+1, x)+(l+m-1)*(l+m)*legendrePolyA(l-1,m-1,x))

	cache = legendrePolyA_cache.get(l, None)
	if(cache == None): # w/o l
		legendrePolyA_cache[l] = {}
		legendrePolyA_cache[l][m] = {}
		legendrePolyA_cache[l][m][x] = ((1-x**2)**(1/2))*(-1)/(2*m)*(legendrePolyA(l-1, m+1, x)+(l+m-1)*(l+m)*legendrePolyA(l-1,m-1,x))
	else: # w l
		cache = legendrePolyA_cache[l].get(m, None)
		if(cache == None): # w/o m
			legendrePolyA_cache[l][m] = {}
			legendrePolyA_cache[l][m][x] = ((1-x**2)**(1/2))*(-1)/(2*m)*(legendrePolyA(l-1, m+1, x)+(l+m-1)*(l+m)*legendrePolyA(l-1,m-1,x))
		else:
			cache = legendrePolyA_cache[l].get(x, None)
			if(cache == None): # w/o x
				legendrePolyA_cache[l][m][x] = ((1-x**2)**(1/2))*(-1)/(2*m)*(legendrePolyA(l-1, m+1, x)+(l+m-1)*(l+m)*legendrePolyA(l-1,m-1,x))

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
# FIXME: 45 items is also an arbitrary precision limit.
def cos(theta):
	result = 1
	for i in range(1, 45):
		result += (-1)**i * theta**(i*2) / factorial(i*2)
	return result

# Utility: Exponential Function
# Uses Taylor Series Expansion
# Depends on:
#    factorial
# FIXME: 32 items an arbitrary precision limit
def exp(x):
	result = 1
	for i in range(1, 32):
		result += x**i / factorial(i)
	return result

# SphericalHarmonics(l, m, theta, phi)
# Computes spherical harmonics Ylm(theta, phi) for given parameters.
# Depends on the following functions:
#    legendrePolyA
#    cos (cosine function)
#    exp (exponential function)
#
def SphericalHarmonics(l, m, theta, phi):
	return ((2*l+1)*factorial(l-m)/factorial(l+m)/(4*pi))**(1/2) * legendrePolyA(l, m, cos(theta)) * exp(1j * m * phi)


# print(legendrePolyNA(5000, pi/1000))
print(SphericalHarmonics(100, 1, pi/1000, pi/6))