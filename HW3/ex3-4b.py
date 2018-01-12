####################################################################
# Computational Physics, 2017-18 Sem1
# HW-3 Ex-3-4b (Drawing Plots)
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
xMin = -600
xMax =  600
dx   =  0.1

# t range
tMin = 0
tMax = 221 # 2Npi/w
dt   = 0.1

# Configure Electric Field time-dependency
E   = lambda t: 0.05338027295 * (sin(0.01423854766 * t))**2 * sin(0.1139083813 * t)

###################################################################
# Matrix Generators 
# Most Arithmetics have been pre-computed for faster initialization.
NX   = int((xMax - xMin)/dx)
NT   = int((tMax - tMin)/dt)

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

#################  PREGENERATED EX3-4A RESULT  ###################
# !!! CHECK DOCUMENTATION TO KNOW THE PARAMETERS FOR GENERATION !!!
# !!! YOU CAN CALCULATE THE SAME RESULTS USING EX3-4.PY !!!
Ts = [tMin + dt * i for i in range(NT)]
asT = [complex(line.strip()) for line in open("ex3-4-asT.dat", 'r')]

# # For integration DEBUG ONLY -- use ex3-3b.py
# # 1/sqrt(2pi) = 0.398942280401433
Aa = lambda t0, w: sum([e**(-1j * w * Ts[i]) * e**(-(Ts[i] - t0)**2/(2*15*15)) * asT[i] for i in range(NT)]) * dt

from math import log10
markerws = [0.1 * (1+i) for i in range(210)]
realws  = [marker * 0.1139083813 for marker in markerws]
t0s = [0.2 * (1+j) for j in range(1100)]

# Aas = [log10(abs(Aa(w))**2) for w in realws]
Aas = [[log10(abs(Aa(t0, w))**2) for w in realws] for t0 in t0s]

import matplotlib.pyplot as plt
fig, ax = plt.subplots()
# contourf(X, Y, Z)
# len(Y) must be # of rows in Z, so rows of Z are grouped by values so
# Z = [[value(x1, y1), value(x2, y1), ...], [value(x1, y2), ...]]
ax.contourf(markerws, ts, Aas, cmap="coolwarm")
plt.show()