# $ RELEASE $
# $ 201711060205Z $ rel01
# $ Signed-Off-By: Haipeng Lin <jimmie.lin@gmail.com>
####################################################################
# Computational Physics, 2017-18 Sem1
# HW-1 Ex-2
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
# Now Playing: I know - Tom Odell
####################################################################

# Solves a linear differential equation with an initial condition
# dy/dx = f(x, y)  a <= x <= b
# y(a) = y_0

# Euler Forward Method
# eulerForwardSolve(\x, y -> f(x,y), a, b, y0, dt = None)
# If timestep dt is none, then a default value will be used given N = 1000 partitions
# Returns the values for all the partitioned f(x_i) = y_i
#
# FIXME: 1000 partitions is completely arbitrary
# although not relevant in this exercise as dt is given in the problem.
def eulerForwardSolve(f, a, b, y0, dt = None):
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

    return ys

# Actual Code for exercise 2
# Solve
#   dN(t)        N(t)
#  ------- = - -------
#    dt          tau       (tau = 1s)
#
print("dt = 0.4s")
r1 = (eulerForwardSolve(lambda x, y: (-1)*y/1, 0, 5, 100, 0.4))
for i in range(1, len(r1) + 1):
    print(r1[i-1])

print("dt = 0.2s")
r2 = (eulerForwardSolve(lambda x, y: (-1)*y/1, 0, 5, 100, 0.2))
for i in range(1, len(r2) + 1):
    print(r2[i-1])

print("dt = 0.1s")
r3 = (eulerForwardSolve(lambda x, y: (-1)*y/1, 0, 5, 100, 0.1))
for i in range(1, len(r3) + 1):
    print(r3[i-1])

print("dt = 0.05s")
r4 = (eulerForwardSolve(lambda x, y: (-1)*y/1, 0, 5, 100, 0.05))
for i in range(1, len(r4) + 1):
    print(r4[i-1])
