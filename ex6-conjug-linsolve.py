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

####################################################################
# The Lightweight Matrix Library, Python 3+
# (c) 2017 Haipeng Lin <linhaipeng@pku.edu.cn>
# Version 1710.01
#
# edu.jimmielin.pylitematrix
#
# This library is fully written by myself to provide for some basic
# matrix operation functions for homework writing.
# It is included inline for readability.
#
# Loosely based on the Haskell Matrix Functions written by myself,
# used for homework assignments in the "Numerical Methods B" course
# in the 2016/17 Semester 1.
# The paradigms and crude variable naming may look like Haskell,
# e.g. fmap, constructor functions, excessive use of lambdas and
# recursion, etc...
# <3 Haskell (sadly its not allowed in our Homework Assignments)
####################################################################

## Matrix Generators
# matrix(rows, cols, \i, j -> elem)
# Returns a matrix with data filled according to the constructor lambda.
def matrix(r, c, fn):
    matrix = []
    for i in range(1, r+1):
        row = []
        for j in range(1, c+1):
            row.append(fn(i, j))
        matrix.append(row)
    return matrix

# zeros(i, j)
# Returns a zero filled matrix
def zeros(i, j):     return matrix(i, j, lambda i, j: 0)

# diag(n, list)
# Returns a diagonal matrix with first n elems from list taken to construct
def diag(n, list):   return matrix(n, n, lambda i, j: list[i] if i == j else 0)

## Matrix Manipulation Functions
# mmap(matrix, fn)
# equiv. fmap fn matrix
def mmap(matrix, fn):
    for i, row in enumerate(matrix):
        for j, val in enumerate(row):
            val = fn(val)
    return matrix

# dotproduct_(list1, list2)
# Multiplies two *lists* of the *same* length
def dotProduct_(l1, l2):
    if(len(l1) != len(l2)):
        raise ValueError("Attempted to multiply two different-sized lists", l1, l2)
    x = 0
    for i, n in enumerate(l1):
        x += l1[i] * l2[i]
    return x

# Multiplies two *matrices* representing row/col vectors
def dotProduct(m1, m2):
    return dotProduct_(sum(m1, []), sum(m2, [])) # good taste

# nrows(matrix)
# Gets number of rows
def nrows(matrix):   return len(matrix)

# ncols(matrix)
# Gets number of cols
def ncols(matrix):   return (0 if len(matrix) == 0 else len(matrix[0]))

# getRow(matrix, r) O(1)
def getRow(matrix, r):
    if(len(matrix) < r):
        raise ValueError("Out-of-Bounds Matrix Row Access", r, matrix)
    return matrix[r-1]

# getCol(matrix, c) O(r)
def getCol(matrix, c):
    col = []
    if(len(matrix) == 0 and c > 0):
        raise ValueError("Out-of-Bounds Matrix Column Access", c, matrix)
    elif(len(matrix[0]) < c):
        raise ValueError("Out-of-Bounds Matrix Column Access", c, matrix)

    for i, row in enumerate(matrix):
        col.append(row[c-1])
    return col

# elem(matrix, i, j)
# Get [i, j] in matrix
def elem(matrix, i, j):
    if(len(matrix) < i or len(matrix[i-1]) < j):
        raise ValueError("Out-of-Bounds Matrix Access", i, j, matrix)
    return matrix[i-1][j-1]

## Matrix Computational Functions
# mult(matrix1, matrix2)
# multStd__ a b = matrix r c $ \(i,j) -> dotProduct (V.unsafeIndex avs $ i - 1) (V.unsafeIndex bvs $ j - 1)
# multStd_ a@(M n m _ _ _ _) b@(M _ m' _ _ _ _) = matrix n m' $ \(i,j) -> sum [ a !. (i,k) * b !. (k,j) | k <- [1 .. m] ]
# multStd is generally slower than multStd2 but if you need debugging, multStd_ is more by definition
def mult(m1, m2):
    return matrix(nrows(m1), ncols(m2), lambda i, j: dotProduct_(getRow(m1, i), getCol(m2, j)))

####################################################################
# / end Lightweight Matrix Library
####################################################################

def matrixA(n):
    matrix = []
    for i in range(1,n+1):
        row = []
        for j in range(1,n+1):
            row.append(-1 if (abs(j-i) == 1) else (3 if j == i else (0.5 if j == (n+1-i) else 0)))
        matrix.append(row)
    return matrix

def vectorB(n):     # assuming n is even as implied by the definition (2.5, 1.5...1.5, 1.0, 1.5...2.5)T
    vector = []

    return vector