# $ RELEASE $
# $ 201711060212Z $ rel01
# $ Signed-Off-By: Haipeng Lin <jimmie.lin@gmail.com>
####################################################################
# Computational Physics, 2017-18 Sem1
# HW-1 Ex-9
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
# Now Playing: 野孩子 - 杨千嬅
####################################################################

# Euler's Number Constant
e = 2.718281828459045

# This is a simple, Python-dependent implementation of the exp(x)
# If you want to see a Taylor expansion based implementation
# (which is usually more accurate than the Pythonic one)
# Check out Ex-4's Source Code
def exp(x):
    return e**x

# Trapezoid Formula for Integration
# Num integTrapezoid(f lambda x: y, a, b, h)
def integTrapezoid(f, a, b, n):
    result = 0
    h  = (b - a)/n

    # Sum over k = 1, ... n-1 for f(x_k)
    result = h/2 * (f(a) + f(b) + 2*sum([f(a+i*h) for i in range(1,n)]))

    return result

# Ex-9 Specific Code
f = lambda r: (-1)/2187 * r**4 * exp((-1) * 2 * r / 3) * (4/81*r**2 - 16/27*r + 16/9)
print(integTrapezoid(f, 0, 60, 1000))