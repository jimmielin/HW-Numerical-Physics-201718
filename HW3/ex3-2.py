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
xMin = -2000
xMax =  2000
dx   =     0.1

# t range
tMin = 0
tMax = 200 # 2Npi/w
dt   = 0.05  

# Configure Electric Field time-dependency
E   = lambda t: 0.016880324465 * (sin(0.025 * t))**2 * sin(t)

######### NO USER CONFIGURABLE OPTIONS BEYOND THIS POINT ##########

####################################################################
# UTILITY FUNCTIONS
####################################################################
pi = 3.141592653589793
e  = 2.718281828459045

# Float|Complex sin(Float|Complex theta)
def sin(theta):
    r = 0.5*(e**(1j*theta) - e**(-1j*theta))/1j
    return (r.real if not isinstance(theta, complex) else r)

# Float|Complex cos(Float|Complex theta)
def cos(theta):
    r = 0.5*(e**(1j*theta) + e**(-1j*theta))
    return (r.real if not isinstance(theta, complex) else r)

#################  PREGENERATED 1S EIGENVECTOR  ###################
# !!! CHECK DOCUMENTATION TO KNOW THE PARAMETERS FOR GENERATION !!!
# !!! YOU CAN CALCULATE THE SAME RESULTS USING EX3-1.PY !!!
Psi1S = [float(line.strip()) for line in open("ex3-1-eigenvector.dat", 'r')]

###################################################################
# Matrix Generators 
# Most Arithmetics have been pre-computed for faster initialization.
NX   = int((xMax - xMin)/dx)
NT   = int((tMax - tMin)/dt)

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
# Step 2 - (2)
# Calculate:
#   a) \Psi_1s(x)'s time evolution to \Psi(x, t)
#   b) P(k) = |<\Psi_k | \Psi_f>|^2 \Psi_f=\Psi(t_f)-p\Psi_1s, p=<\Psi_1s|\Psi(tf)>
#      k \in [-1.5, 1.5]
# Tridiag Matrix is stored as a, b, c

# Generate Matrices and base data...
dxISq = 1/dx**2

# Starting condition
# Psi_1s is vector
# Time evolved vectors are stored in Psis in 2-D tuple (t, [...])
P1ss = [(0, Psi1S)]

# Evolve through time. Advances through timestep x, given P1s (from P1ss) in mapping (t, Psi) 2d-tuple
# This *reconstructs* the hamiltonian because each different component must be previously
# fully-evolved before we start the next one.
def nextPsi(P1s):
    Psi = P1s[1][:]
    t   = P1s[0]
    tp  = P1s[0] + dt/2

    # Psi2
    # Solve tridiagonal matrix for Psi2, where
    # A * Psi2 = B * Psi

    # Generators for A
    dxISq = 1/dx**2
    Abs = [(dxISq - 1/(1 + (xMin + i * dx)**2)**(1/2) + E(t + dt/2) * (xMin + i * dx)) * (0.5 * 1j * dt) + 1 for i in range(NX)]

    Ac = lambda i: -0.5 * dxISq * (0.5 * 1j * dt)
    Aa = lambda i: Abs[i-1]
    Ab = lambda i: -0.5 * dxISq * (0.5 * 1j * dt)

    # Generators for B
    Bbs = [(dxISq - 1/(1 + (xMin + i * dx)**2)**(1/2) + E(t + dt/2) * (xMin + i * dx)) * (-0.5 * 1j * dt) + 1 for i in range(NX)]

    Bc = lambda i: -0.5 * dxISq * (-0.5 * 1j * dt)
    Ba = lambda i: Bbs[i-1]
    Bb = lambda i: -0.5 * dxISq * (-0.5 * 1j * dt)

    # Multiply to get right-vector
    target = triMultVector(Ba, Bb, Bc, Psi)

    # Solve
    Psi2 = triSolve(Aa, Ab, Ac, target)

    # DIAGNOSTICS: check solution stability
    # targett = triMultVector(Aa, Ab, Ac, Psi2)
    # offset = max([abs(targett[i] - target[i]) for i in range(len(targett))])
    # print("* diag nextPsi: chk =", offset, flush=True)

    return (P1s[0] + dt, vecNormalize(Psi2))

# For debugging only: Echoes some debug information for given tuple (t, Psi)
def debugP1s(P1s):
    Psi = P1s[1]
    print("t =", P1s[0], "Psi excerpt:", Psi[0:3], "...", Psi[(-4):])

currT = 0
P1sT = [abs(vecDotProduct(P1ss[0][1], P1ss[-1][1]))**2] # 布居数随时间的演化
debugP1s(P1ss[0])
print(P1ss[-1][0], P1sT[-1])
for i in range(1, NT):
    P1ss.append(nextPsi(P1ss[-1]))
    # debugP1s(P1ss[-1])
    # Dot product to get the result...
    P1sT.append(abs(vecDotProduct(P1ss[0][1], P1ss[-1][1]))**2)
    print(P1ss[-1][0], P1sT[-1])

import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ts = [0 + i * dt for i in range(NT)]
ax.plot(ts, P1sT)
plt.show()

print(P1ss[-1])