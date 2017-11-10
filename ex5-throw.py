####################################################################
# Computational Physics, 2017-18 Sem1
# HW-1 Ex-5
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
# Now Playing: My Heart Will Go On - Celine Dion
#              Not Another Song About Love - Hollywood Ending
####################################################################

# The below is mostly adapted from the code in Exercise 2
# Solves a linear differential equation with an initial condition
# dy/dx = f(x, y)  a <= x <= b
# y(a) = y_0

# Euler Forward Method, for XY coupled situation
# eulerForwardSolve(\t, yx, yy -> f(t,yx), \t, yy, yx -> f(t,yy), a, b, x0, y0, dt = None)
# If timestep dt is none, then a default value will be used given N = 1000 partitions
# Returns the values for all the partitioned f(t_i) in a tuple (t_i, x_i, y_i)
#
# FIXME: default 1000 partitions is completely arbitrary
def eulerForwardSolveXY(fx, fy, a, b, x0, y0, dt = None):
	# Euler forward method
    k = 0
    if(dt != None):
        h = dt    # timestep given
    else:
        h = (b-a)/1000
    xn = a
    ys = [(a,x0,y0)]

    while((xn) <= b): # do a little nudging because python's double is a little too sensitive
        xn += h
        ys.append((xn, ys[-1][1] + h * fx(xn, ys[-1][1], ys[-1][2]), ys[-1][2] + h * fy(xn, ys[-1][2], ys[-1][1])))

    return ys


# Euler Forward Method
# eulerForwardSolve(\x, y -> f(x,y), a, b, y0, dt = None)
# If timestep dt is none, then a default value will be used given N = 1000 partitions
# Returns the values for all the partitioned f(x_i) = y_i
#
# FIXME: 1000 partitions is completely arbitrary
def eulerForwardSolve(f, a, b, y0, dt = None, terminateCondition = None):
    k = 0
    if(dt != None):
        h = dt    # timestep given
    else:
        h = (b-a)/1000
    xn = a
    ys = [(a,y0)]

    while((xn) <= b): # do a little nudging because python's double is a little too sensitive
        xn += h
        ys.append((xn, ys[-1][1] + h * f(xn, ys[-1][1])))
        if(terminateCondition != None): # then is a lambda
        	if(terminateCondition(xn, ys[-1][1]) == True):
        		break

    return ys

# The below is (mostly) copy-pasted from the code in Exercise 4
# Trigonometry Functions
pi = 3.141592653589793

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
# Uses Taylor Series Approximation for sine
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

def deg2Rad(deg):
	return deg * 0.017453292519943295

## Exercise 5 Code
# Configure parameters below:
# degTheta = 30      # deg
v0 = 700           # m/s
hasFriction = True
fFactor = 4e-5     # M^(-1) as friction coefficient
g = 9.8            # m/s^(-2) as gravity constant
dt = 0.005          # dt as timestep
te = 10000         # Arbitrarily defined maximum time_end (not actually used, as lambda t, y: y <= 0 termination rule will)
#                  # control the result from breaking the floor, heh

if(hasFriction == False):
	fFactor = 0
	print("** WARNING: Friction set to 0 **")

degThetas = [30, 40, 50]
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
for i in range(len(degThetas)):
	degTheta = degThetas[i]
	print("Calculating for theta = ", degTheta, "deg")

	# Solve for dvy/dt, dvx/dt first
	vs = (eulerForwardSolveXY(lambda t, vx, vy: (-1) * fFactor * (vx**2 + vy**2) ** (1/2) * vx, lambda t, vy, vx: (-1) * g - fFactor * (vx**2 + vy**2) ** (1/2) * vy, 0, te, v0 * cos(deg2Rad(degTheta)), v0 * sin(deg2Rad(degTheta)), dt))

	# Then simply solve for x and y...
	# Translate vs to an appropriate data structure (dictionary)
	vxs = {}
	vys = {}
	for i in range(len(vs)):
		vxs[vs[i][0]] = vs[i][1]
		vys[vs[i][0]] = vs[i][2]

	xs_c = (eulerForwardSolve(lambda t, x: vxs[t], 0, te, 0, dt))
	ys_c = (eulerForwardSolve(lambda t, y: vys[t], 0, te, 0, dt, lambda t, y: y <= 0))

	xs = []
	ys = []

	for i in range(len(ys_c)):
		xs.append(xs_c[i][1])
		ys.append(ys_c[i][1])

	# Slice xs to fit ys...
	xs = xs[0:len(ys)]

	# Compute parameters as required by the exercise.
	print("* Fall time t =", xs_c[len(ys_c)-1][0], "s")
	print("* Velocity vx =", vs[len(ys_c)-1][1], "m/s, vy =", vs[len(ys_c)-1][2], "m/s")
	print("* Speed v =", (vs[len(ys_c)-1][1]**2 + vs[len(ys_c)-1][2]**2)**(1/2), "m/s")
	print("* Distance xf =", xs_c[len(ys_c)-1][1], "m")
	print("\n")

	ax.plot(xs, ys, label=degTheta.__str__())

hasFriction = False
fFactor = 0
for i in range(len(degThetas)):
    degTheta = degThetas[i]
    print("Calculating for theta = ", degTheta, "deg")

    # Solve for dvy/dt, dvx/dt first
    vs = (eulerForwardSolveXY(lambda t, vx, vy: (-1) * fFactor * (vx**2 + vy**2) ** (1/2) * vx, lambda t, vy, vx: (-1) * g - fFactor * (vx**2 + vy**2) ** (1/2) * vy, 0, te, v0 * cos(deg2Rad(degTheta)), v0 * sin(deg2Rad(degTheta)), dt))

    # Then simply solve for x and y...
    # Translate vs to an appropriate data structure (dictionary)
    vxs = {}
    vys = {}
    for i in range(len(vs)):
        vxs[vs[i][0]] = vs[i][1]
        vys[vs[i][0]] = vs[i][2]

    xs_c = (eulerForwardSolve(lambda t, x: vxs[t], 0, te, 0, dt))
    ys_c = (eulerForwardSolve(lambda t, y: vys[t], 0, te, 0, dt, lambda t, y: y <= 0))

    xs = []
    ys = []

    for i in range(len(ys_c)):
        xs.append(xs_c[i][1])
        ys.append(ys_c[i][1])

    # Slice xs to fit ys...
    xs = xs[0:len(ys)]

    # Compute parameters as required by the exercise.
    print("* Fall time t =", xs_c[len(ys_c)-1][0], "s")
    print("* Velocity vx =", vs[len(ys_c)-1][1], "m/s, vy =", vs[len(ys_c)-1][2], "m/s")
    print("* Speed v =", (vs[len(ys_c)-1][1]**2 + vs[len(ys_c)-1][2]**2)**(1/2), "m/s")
    print("* Distance xf =", xs_c[len(ys_c)-1][1], "m")
    print("\n")

    ax.plot(xs, ys, label=degTheta.__str__() + " (Friction off)")

plt.legend(loc='upper left')
plt.show()