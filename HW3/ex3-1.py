####################################################################
# Computational Physics, 2017-18 Sem1
# HW-3 Ex-3
#
# (c) 2017-2018 Haipeng Lin <linhaipeng@pku.edu.cn>
# All Rights Reserved.
#
# This program is written as Homework for the Computational Physics
# Course in Peking University, School of Physics.
#
# TO THE MAXIMUM EXTENT PERMITTED BY LAW,
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR "AS IS", WITHOUT ANY EXPRESS
# OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
# IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, 
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
# IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF 
# THE POSSIBILITY OF SUCH DAMAGE.
#
# Copying is strictly prohibited and its usage is solely restricted
# to uses permitted by the author and homework grading.
#
# This program is Python 3 (3.6.2 on MSYS_NT)
# compatible, and does not support Python 2.
#
# Now Playing: Where We Belong - Thomas Fiss
####################################################################

####################################################################
# CONFIGURABLE PARAMETERS
# For HW-3 Ex-3
####################################################################

# x range
xM   = 2000
xMin = -xM
xMax =  xM
dx   = 0.1

# Ansatz for eigenvalue of 1s
offset = 0.67

######### NO USER CONFIGURABLE OPTIONS BEYOND THIS POINT ##########

###################################################################
# Matrix Generators 
# Most Arithmetics have been pre-computed for faster initialization.
NX   = int((xMax - xMin)/dx)

###################################################################
# Tridiagonal Matrix Computational Functions
# Provide lambdas a, b, c defining tridiagonal matrices as
# (a1 b1                          )
# (c1 a2 b2                       )
# (   c2 a3 b3                    )
# (      .. .. ..                 )
# (                           bn-1)
# (                       cn-1  an)

# Multiply by vector
# List triMultVector(lambda a, lambda b, lambda c, List vector)
def triMultVector(a, b, c, vector):
    r = len(vector)
    flat = [a(1) * vector[0] + b(1) * vector[1]]
    for rr in range(1, r - 1):
        flat.append(c(rr) * vector[rr-1] + a(rr+1) * vector[rr] + b(rr+1) * vector[rr+1])
    flat.append(c(r-1) * vector[r-2] + a(r) * vector[r-1])

    return flat

# Solve tridiagonal matrix (adapted from HW-1 Ex-8, but uses lambdas)
# Ax = y
# Using tridiagonal matrix algorithm ("Thomas algorithm")
def triSolve(a, b, c, y):
    n = len(y)

    # The arrays below have been padded so the indeces are mathematically sound.
    # One doesn't have to do this, but it makes for more coding convenience.
    A = [0, 0] # lower diagonal A, i=2,3..n
    B = [0] # diagonal B, i=1,2,..n
    C = [0] # upper diagonal C, i=1,2,..n-1
    D = y[:] # right-b vector d1, d2, ... n
    D.insert(0, 0) # left-padd, * minor performance concern

    # The solution array is *NOT* padded, as it'll be used for
    # fromList_ by the Matrix prototype
    # It's even reversed as part of back-substitution
    # so watch out for the index.
    X = [None] * n

    # .elem is now O(1) so you can use that for safety,
    # but I'll tap into .internal represenation for speed and less overhead
    # since we already checked the indexes above
    for i in range(1, n+1):
        if(i != 1): # make sure i=1 case does not run for A
            A.append(c(i))
        if(i != n): # make sure i=n case does not run for C
            C.append(b(i))
        B.append(a(i))

    # Now perform the Thomas Algorithm Forward Sweep
    C[1] = C[1] / B[1]
    D[1] = D[1] / B[1]
    for i in range(2, n+1):
        if(i != n): # make sure i=n case does not run for C=..n-1
            C[i] = C[i] / (B[i] - A[i] * C[i-1])
        D[i] = (D[i] - A[i] * D[i-1]) / (B[i] - A[i] * C[i-1])

    # Now do a back-subsitution...
    X[0] = D[n]
    for i in range(2, n+1):
        j = n+1-i # back indexes
        X[i-1] = D[j] - C[j] * X[i-2]

    X.reverse()

    return X

# "Dot Product" two Vector-Lists
def vecDotProduct(a, b):
    if(len(a) != len(b)):
        raise TypeError("vecDotProduct: cannot dot product two different dimensional lists")

    return sum([a[i].conjugate() * b[i] for i in range(len(a))]) * dx

# Convert Vector to Unity
def vecNormalize(vector):
    Norm = vecDotProduct(vector, vector) ** (1/2)
    vecNormd   = [vectorx/Norm for vectorx in vector] # this is the normalized eigenvector

    return vecNormd

###################################################################
# Step 1 - Compute Hydrogen Ground State
# Tridiag Matrix is stored as a, b, c
# Use A - (-0.5I) as approx is -0.5
dxISq = 1/dx**2
diags = [(dxISq - (1 + (xMin + i * dx)**2)**(-1/2) + offset) for i in range(NX)]

c = lambda i: -0.5 * dxISq
a = lambda i: diags[i-1]
b = lambda i: -0.5 * dxISq

# Testing code for triMultVector & triSolve
# testa = lambda i: 2
# testb = lambda i: -1
# testc = lambda i: -1
# vectest = [1 for i in range(5)]
# print(vectest)
# x = triSolve(testa, testb, testc, vectest)
# print(x)
# v = triMultVector(testa, testb, testc, x)
# print(v)

# vectest = [1 for i in range(5)]
# print(vectest)
# x = triSolve(testa, testb, testc, vectest)
# print(x)
# v = triMultVector(testa, testb, testc, x)
# print(v)
# z = [0.5 + i for i in range(1, 500)]
# print(z)
# x = triSolve(lambda i: 2.38, lambda i: -1, lambda i: 0.392, z)
# zz = triMultVector(lambda i: 2.38, lambda i: -1, lambda i: 0.392, x)
# print(zz)

# exit()

# Solve using inverse power method:
# Ay_k = z_(k-1)
#  u_k = ||y_k|| infty
#  z_k = y_k / u_k
# The first step is basically solving for y_k, in relation to a
# tridiagonal matrix.
# It can be solved using the thomas transform method ("追赶法")
z = [1 for i in range(NX)] # Just a starting vector, bad guess...
for i in range(1, 100):
    y = triSolve(a, b, c, z)
    yInftyMod = (max(-min(y), max(y))) # use abs. value

    ## DIAGNOSTICS: Perform a consistency check (this is slow, remove in production)
    # zz = triMultVector(a, b, c, y)
    # df = [zz[i] - z[i] for i in range(len(z))]
    # print("* diff:", max(df))
    # print("*", end='', flush=True)
    print(1/yInftyMod - offset, flush=True)

    z = [yi/yInftyMod for yi in y]

print("")
print("Lambda =", 1/yInftyMod)
print("E0 =", 1/yInftyMod - offset)

# Get the state (z) which is an eigenvector and compute its euclidean length -- note normalization
# PSq   = [p**2 for p in z]
# Norm  = ((sum(PSq)) * dx) * (1/2)
# Psi   = [zx/Norm for zx in z] # this is the normalized eigenvector
Psi = vecNormalize(z)
print(vecDotProduct(z, z))
print(vecDotProduct(Psi, Psi))

PsiSq = [PsiX**2 for PsiX in Psi]
xs    = [xMin + i * dx for i in range(NX)]

# Draw PsiSq across -x ~ x axes
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.plot(xs, PsiSq)
plt.show()

print(Psi)