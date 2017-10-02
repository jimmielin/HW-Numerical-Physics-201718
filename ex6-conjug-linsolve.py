####################################################################
# Computational Physics, 2017-18 Sem1
# HW-1 Ex-8
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
# Now Playing: Forgettable - Project 46
####################################################################

import math
import time

# Include my own lightweight matrix library
from pylitematrix.pylitematrix import Matrix

def matrixA(n): 
    return Matrix(n, n, lambda i, j: -1 if (abs(j-i) == 1) else (3 if j == i else (0.5 if j == (n+1-i) else 0)))

def vectorB(n): 
    if(n % 2):
        return Matrix(n, 1, lambda m, _: 2.5 if(m == 1 or m == n) else (1.0 if (m == math.ceil(n/2)) else 1.5))
    else:
        return Matrix(n, 1, lambda m, _: 2.5 if(m == 1 or m == n) else (1.0 if (m == n/2 or m == n/2+1) else 1.5))

# The Conjugate Gradient Method, w/o Improved Algorithm
def conjGradient(A, b, x0):
    rc = b - A * x0
    if(rc.max() == 0):
        return x0
    pc = rc
    k  = 0
    x  = x0
    while(k < 1000):
        alpha_factor = (rc.transpose() * A * pc).elem(1, 1)
        if(alpha_factor == 0):
            break
        alpha = (rc.transpose() * pc).map(lambda x: x/alpha_factor).elem(1, 1)
        x = x + pc.map(lambda x: x * alpha)
        beta_factor  = (pc.transpose() * A * pc).elem(1, 1) * -1
        rc = rc - (A * pc).map(lambda x: x * alpha)
        if(beta_factor == 0):
            break
        beta = (rc.transpose() * A * pc).map(lambda x: x/beta_factor).elem(1, 1)
        pc = rc + pc.map(lambda x: x * beta)
        k = k + 1
    print("debug from conjGradient: did ", k, "iterations")
    return x

# The Conjugate Gradient Method, w/ Improved Algorithm
# Generally is 50% faster than conjGradient
def conjGradient2(A, b, x0):
    rc = b - A * x0
    if(rc.max() == 0):
        return x0
    pc = rc
    k  = 0
    x  = x0
    while(k < 1000):
        # print("*", end="", flush=True) # iterator diagnostic
        alpha_factor = (pc.transpose() * A * pc).elem(1, 1)
        if(alpha_factor == 0):
            break
        alpha = (rc.transpose() * rc).map(lambda x: x/alpha_factor).elem(1, 1)
        x = x + pc.map(lambda x: x * alpha)
        beta_factor  = (rc.transpose() * rc).elem(1, 1)
        rc = rc - (A * pc).map(lambda x: x * alpha)
        if(beta_factor == 0):
            break
        beta = (rc.transpose() * rc).map(lambda x: x/beta_factor).elem(1, 1)
        pc = rc + pc.map(lambda x: x * beta)
        k = k + 1
    print("debug from conjGradient2: did ", k, "iterations")
    return x

# Actual Code for Exercise 6
def ex6(n):
    print(conjGradient(matrixA(n), vectorB(n), Matrix(n, 1, lambda i, j: 0)))

def ex6_2(n):
    print(conjGradient2(matrixA(n), vectorB(n), Matrix(n, 1, lambda i, j: 0)))

def ex6_2_noprint(n):
    (conjGradient2(matrixA(n), vectorB(n), Matrix(n, 1, lambda i, j: 0)))

timings = []
for i in range(2, 1000):
    start_time = time.time()
    ex6_2_noprint(i) # don't print
    elapsed_time = time.time() - start_time
    timings.append(elapsed_time)
    print("Execution time is i=", i, " t=", elapsed_time, " seconds")

print(timings)

#start_time = time.time()
#ex6(100)
#elapsed_time = time.time() - start_time
#print("Execution time is ", elapsed_time, " seconds.")

#start_time = time.time()
#ex6_2(10000)
#elapsed_time = time.time() - start_time
#print("Execution time is ", elapsed_time, " seconds.")