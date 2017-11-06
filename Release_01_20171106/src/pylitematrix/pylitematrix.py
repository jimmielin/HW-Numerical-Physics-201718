####################################################################
# The Lightweight Matrix Library, Python 3+
# (c) 2017 Haipeng Lin <linhaipeng@pku.edu.cn>
# Version 1710.02
#
# edu.jimmielin.pylitematrix
#
# This library is fully written by myself to provide for some basic
# matrix operation functions for homework writing.
#
# Loosely based on the Haskell Matrix Functions written by myself,
# used for homework assignments in the "Numerical Methods B" course
# in the 2016/17 Semester 1.
# The paradigms and crude variable naming may look like Haskell,
# e.g. fmap, constructor functions, excessive use of lambdas and
# recursion, etc...
# <3 Haskell (sadly its not allowed in our Homework Assignments)
####################################################################

import math

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
    def diag_(self, n, list):   return Matrix(n, n, lambda i, j: list[i-1] if i == j else 0)

    # fromList(rows, cols, list)
    # From a list with at least rows*cols elements generate a matrix left-to-right, top-to-down
    # fromList n m = M n m 0 0 m . V.fromListN(n*m)
    def fromList_(self, r, c, list):
        if(len(list) < r*c):
            raise ValueError("List too short for fromList Matrix Construction", r*c, len(list))
        return Matrix(r, c, lambda i, j: list[(i-1)*c+j-1])

    ## Matrix Manipulation Functions
    # map(fn)
    # equiv. fmap fn matrix, return a copy (for in-place changes, use mapInplace, which is map! ruby)
    def mapInplace(self, fn):
        for i, row in enumerate(self.internal):
            for j, val in enumerate(row):
                self.internal[i][j] = fn(val) # because we are not by reference
        return self

    def map(self, fn):
        mp = Matrix(1, 1, lambda i, j: 0)
        mp.internal = [xs[:] for xs in self.internal]
        # deepcopy because here it is a reference now by default, stupid languages with side effects
        mp.mapInplace(fn)
        return mp

    # toList(matrix)
    # Get elements of matrix stored in a list.
    def toList(self):  return sum(self.internal, [])

    # dotproduct_(list1, list2)
    # Multiplies two *lists* of the *same* length
    def dotProduct__(self, l1, l2):
        if(len(l1) != len(l2)):
            raise ValueError("Attempted to multiply two different-sized lists", l1, l2)
        x = 0
        for i, n in enumerate(l1):
            x += l1[i] * l2[i]
        return x

    # Multiplies two *matrices* representing row/col vectors (Prototype)
    def dotProduct(self, m1, m2):
        return self.dotProduct__(sum(m1, []), sum(m2, [])) # good taste

    # nrows(matrix)
    # Gets number of rows (Prototype)
    def nrows_(self, m):   return len(m)

    # ncols(matrix)
    # Gets number of cols (Prototype)
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

    # subMatrix(s_r, e_r, s_c, e_c) (starting/ending row/column)
    # extract a submatrix given row and column limits, e.g.
    # subMatrix(1, 2, 2, 3) (1 2 3) = (2 3)
    #                       (4 5 6) = (5 6)
    #                       (7 8 9)
    # FIXME: this function is currently unsafe
    def subMatrix(self, sr, er, sc, ec):
        # M nrows ncols rowOffset colOffset vcols mvect
        # M (r2-r1+1) (c2-c1+1) (ro+r1-1) (co+c1-1) w v
        return Matrix(er-sr+1, ec-sc+1, lambda i, j: self.internal[sr+i-2][sc+j-2])

    # splitBlocks(self, r, c) = (TL, TR, BL, BR) 4-tuple Matrix
    # makes a block-partition matrix using given element (r,c) as reference,
    # with the reference element staying in the bottom-right corner of the first split.
    #                    (             )   (      |      )
    #                    (             )   ( ...  | ...  )
    #                    (    x        )   (    x |      )
    # splitBlocks(r,c) = (             ) = (-------------) , where x = a(r,c)
    #                    (             )   (      |      )
    #                    (             )   ( ...  | ...  )
    #                    (             )   (      |      )
    # Note that some blocks can end up empty. We use the following notation for these blocks:
    # ( TL | TR )
    # (---------)
    # ( BL | BR )
    def splitBlocks(self, r, c):
        return (self.subMatrix(1, r, 1, c), self.subMatrix(1, r, c + 1, self.ncols()), self.subMatrix(r + 1, self.nrows(), 1, c), self.subMatrix(r + 1, self.nrows(), c + 1, self.nrows()))

    # elem(i, j)
    # Get [i, j] in matrix
    def elem(self, i, j):
        if(len(self.internal) < i or len(self.internal[i-1]) < j):
            raise ValueError("Out-of-Bounds Matrix Access", i, j, self.internal)
        return self.internal[i-1][j-1]

    # setElem_(i,j,val)
    # Set [i,j] in matrix (SIDE EFFECTS, but almost everything is these days)
    def setElem_(self, i, j, val):
        if(len(self.internal) < i or len(self.internal[i-1]) < j):
            raise ValueError("Out-of-Bounds Matrix Access", i, j, self.internal)
        self.internal[i-1][j-1] = val
        return self.internal[i-1][j-1]

    # setElem(i, j, val)
    # Set [i,j] in matrix without side effects (returns a copy), immutable
    def setElem(self, i, j, val):
        mp = Matrix(1, 1, lambda i, j: 0)
        mp.internal = [xs[:] for xs in self.internal]
        mp.setElem_(i, j, val)
        return mp

    # transpose()
    # Get matrix transposition. O(c*r)
    def transpose(self):    return Matrix(self.ncols(), self.nrows(), lambda i, j: self.elem(j, i))

    ## Matrix Computational Functions
    # multStd_ a@(M n m _ _ _ _) b@(M _ m' _ _ _ _) = 
    #    matrix n m' $ \(i,j) -> sum [ a !. (i,k) * b !. (k,j) | k <- [1 .. m] ]
    def __mul__(self, m2):
        res = Matrix(self.nrows(), m2.ncols(), lambda i, j: self.dotProduct__(self.getRow(i), m2.getCol(j)))
        return res

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

    ## Matrix Element Functions
    # max() get maximum element from matrix.
    def max(self): return max(self.toList())

    # min() get miminum element from matrix.
    def min(self): return min(self.toList())

# Initialize a Matrix "Prototype"
# Seriously, this is insane - I'm using Javascript paradigms
# in Python.
# But again, the rest is Haskell. You can't get much better than this.
mp = Matrix(1, 1, lambda i, j: 1)

####################################################################
# / end Lightweight Matrix Library
####################################################################