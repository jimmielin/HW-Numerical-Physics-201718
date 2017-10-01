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

class Matrix:
    internal = []

    ## Matrix Generators
    # matrix(rows, cols, \i, j -> elem)
    # Returns a matrix with data filled according to the constructor lambda.
    def __init__(self, r, c, fn):
        matrix = []
        for i in range(1, r+1):
            row = []
            for j in range(1, c+1):
                row.append(fn(i, j))
            matrix.append(row)
        self.internal = matrix

    # The below are PROTOTYPE functions and return new objects
    # zeros(i, j)
    # Returns a zero filled matrix
    def zeros_(self, i, j):     return Matrix(i, j, lambda i, j: 0)

    # diag(n, list)
    # Returns a diagonal matrix with first n elems from list taken to construct
    def diag_(self, n, list):   return Matrix(n, n, lambda i, j: list[i] if i == j else 0)

    # fromList(rows, cols, list)
    # From a list with at least rows*cols elements generate a matrix left-to-right, top-to-down
    # fromList n m = M n m 0 0 m . V.fromListN(n*m)
    def fromList_(self, r, c, list):
        if(len(list) < r*c):
            raise ValueError("List too short for fromList Matrix Construction", r*c, len(list))
        return Matrix(r, c, lambda i, j: list[(i-1)*c+j-1])

    ## Matrix Manipulation Functions
    # map(fn)
    # equiv. fmap fn matrix
    def map(self, fn):
        for i, row in enumerate(self.internal):
            for j, val in enumerate(row):
                self.internal[i][j] = fn(val) # because we are not by reference
        return self

    # toList(matrix)
    # Get elements of matrix stored in a list.
    def toList(self):       return sum(self.internal, [])

    # dotproduct_(list1, list2)
    # Multiplies two *lists* of the *same* length
    def dotProduct_(self, l1, l2):
        if(len(l1) != len(l2)):
            raise ValueError("Attempted to multiply two different-sized lists", l1, l2)
        x = 0
        for i, n in enumerate(l1):
            x += l1[i] * l2[i]
        return x

    # Multiplies two *matrices* representing row/col vectors
    def dotProduct(self, m1, m2):
        return self.dotProduct_(sum(m1, []), sum(m2, [])) # good taste

    # nrows(matrix)
    # Gets number of rows (Universal)
    def nrows_(self, m):   return len(m)

    # ncols(matrix)
    # Gets number of cols (Universal)
    def ncols_(self, m):   return (0 if len(m) == 0 else len(m[0]))

    # nrows(matrix)
    # Gets number of rows
    def nrows(self):   return self.nrows_(self.internal)

    # ncols(matrix)
    # Gets number of cols
    def ncols(self):   return self.ncols_(self.internal)

    # getRow(matrix, r) O(1)
    def getRow(self, r):
        if(len(self.internal) < r):
            raise ValueError("Out-of-Bounds Matrix Row Access", r, self.internal)
        return self.internal[r-1]

    # getCol(matrix, c) O(r)
    def getCol(self, c):
        col = []
        if(len(self.internal) == 0 and c > 0):
            raise ValueError("Out-of-Bounds Matrix Column Access", c, self.internal)
        elif(len(self.internal[0]) < c):
            raise ValueError("Out-of-Bounds Matrix Column Access", c, self.internal)

        for i, row in enumerate(self.internal):
            col.append(row[c-1])
        return col

    # elem(i, j)
    # Get [i, j] in matrix
    def elem(self, i, j):
        if(len(self.internal) < i or len(self.internal[i-1]) < j):
            raise ValueError("Out-of-Bounds Matrix Access", i, j, self.internal)
        return self.internal[i-1][j-1]

    # transpose()
    # Get matrix transposition. O(c*r)
    def transpose(self):    return Matrix(self.ncols(), self.nrows(), lambda i, j: self.elem(j, i))

    ## Matrix Computational Functions
    # multStd_ a@(M n m _ _ _ _) b@(M _ m' _ _ _ _) = 
    #    matrix n m' $ \(i,j) -> sum [ a !. (i,k) * b !. (k,j) | k <- [1 .. m] ]
    def __mul__(self, m2):
        return Matrix(self.nrows(), m2.ncols(), lambda i, j: self.dotProduct_(self.getRow(i), m2.getCol(j)))

    # add
    # O(m*n) rather than naive implementation
    def __add__(self, m2):
        if(self.nrows() != m2.nrows() or self.ncols() != m2.ncols()):
            raise TypeError("Cannot add matrices that are not of the same size!")
        return self.fromList_(self.nrows(), self.ncols(), [a+b for a, b in zip(self.toList(), m2.toList())])

    # invert for unary arithmetic functions
    def __invert__(self): return self.map(lambda x: -1*x)

    # subtract
    def __sub__(self, m2):
        if(self.nrows() != m2.nrows() or self.ncols() != m2.ncols()):
            raise TypeError("Cannot add matrices that are not of the same size!")
        return self.fromList_(self.nrows(), self.ncols(), [a-b for a, b in zip(self.toList(), m2.toList())])

    ## Matrix Pretty-Print
    # adapted from SO#13214809, herein replicated on fair use basis as it
    # is not a crucial component to HW
    def __str__(self):
        s = [[str(n) for n in row] for row in self.internal]
        lens = [max(map(len, col)) for col in zip(*s)]
        form = '\t'.join('{{:{}}}'.format(x) for x in lens)
        table = [form.format(*row) for row in s]
        return ('\n'.join(table))

# Initialize a Matrix "Prototype"
# Seriously, this is insane - I'm using Javascript paradigms
# in Python.
# But again, the rest is Haskell. You can't get much better than this.
mp = Matrix(1, 1, lambda i, j: 1)

####################################################################
# / end Lightweight Matrix Library
####################################################################

def matrixA(n): 
    return Matrix(n, n, lambda i, j: -1 if (abs(j-i) == 1) else (3 if j == i else (0.5 if j == (n+1-i) else 0)))

def vectorB(n): 
    if(n % 2):
        return Matrix(n, 1, lambda m, _: 2.5 if(m == 1 or m == n) else (1.0 if (m == math.ceil(n/2)) else 1.5))
    else:
        return Matrix(n, 1, lambda m, _: 2.5 if(m == 1 or m == n) else (1.0 if (m == n/2 or m == n/2+1) else 1.5))

def conjGradient(A, b, x0):
    return False
