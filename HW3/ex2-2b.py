####################################################################
# Computational Physics, 2017-18 Sem1
# HW-3 Ex-2 B
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
# Now Playing: Seve - Tez Cadey
####################################################################

from pylitematrix.pylitematrix import Matrix

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

####################################################################
# CONFIGURABLE PARAMETERS
# For HW-3 Ex-2
# !!!! THESE CONFIGURABLE OPTIONS ARE NOT ACTUALLY ALL SUPPORTED.
# !!!! JUST LIKE SOME WRF-CHEM MODEL, ONLY A SPECIFIC COMBO WORKS.
# !!!! For homework problems the below parameters given are just fine.
####################################################################

n = 16
h = 0.1
w1 = 2 * sin(11*pi/34) # (16+1)*2 = 34
TU = 2 * pi/w1

# Change L here from 160 ~ + time to observe different behavior, for
# 2) use 160, 4) to 4000 (actually less is already sufficient.)
L = 160 * TU/h
a = 0.25     # alpha
b = 1        # beta

omega = lambda k: 2 * sin(pi*k/(2*n+2))

####################################################################
# MATRIX FUNCTIONS
####################################################################
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
    return x

####################################################################
# MATRIX GENERATORS
####################################################################
Mx = Matrix(n, n, lambda i, j: sin(pi*i*j/(n+1)))
Q = Matrix(n, 1, lambda i, j: 2.828 if i == 11 else 0)
# 4 = 1/sqrt(2/n=16) ~ 

# Padd with "Mathematical Indexes" for convenience...
x = [0] + conjGradient2(Mx, Q, Matrix(n, 1, lambda i, j: 0)).toList()

# INITIAL CONDITION GENERATOR
# DATA STRUCTURE IS AS FOLLOWS -
#   ps[i][t]
#   qs[i][t]
#   dps[i][t] is p' 1st-derivative.
# We store each particle in the 1-st array level because later in the leapfrog method,
# we actually increment by particle (and not by time)
# This requires later flattening of the array.
qs = [[x[i]] if i != 0 and i != (n+1) else [0] for i in range(n+2)]
#ps = [[0][:]] * (n+2)
# nope, python is incredibly stupid
ps = []
for i in range(n+2):
    ps.append([0])
dps = [[qs[i+1][0] + qs[i-1][0] - 2 * qs[i][0] + b * ((qs[i+1][0]-qs[i][0])**3 - (qs[i][0]-qs[i-1][0])**3)] if i != 0 and i != (n+1) else [0] for i in range(n+2)]

for t in range(1, 5):
    qs[0].append(0)
    qs[n+1].append(0)
    for i in range(1, n+1):
        dps[i].append(qs[i+1][t-1] + qs[i-1][t-1] - 2 * qs[i][t-1] + b * ((qs[i+1][t-1]-qs[i][t-1])**3 - (qs[i][t-1]-qs[i-1][t-1])**3))

        # print(len(ps[i]), len(qs[i]), len(dps[i]))
        # List Indexes... and mathematical indexes
        #

        ps[i].append(dps[i][t] * h + ps[i][t-1])
        qs[i].append(ps[i][t] * h + qs[i][t-1])


# print(qs, '\n', '\n', ps, '\n', '\n', dps)
# # print(len(qs), len(ps), len(dps))
# exit()

t = 5
####################################################################
# TIME EVOLUTION OPERATOR
#
# Uses the Leapfrog Method (辛算法) after tips from fellow classmates
# that other forms are numerically unstable.
#
# If t = 4 above you'll get very funny results! It's amazing how
# numerical instability affects the world.
#
# Time is already in TU = 2pi/w1 units.
####################################################################
while t <= L:
    qs[0].append(0)
    qs[n+1].append(0)
    for i in range(1, n+1):
        try:
            ps[i].append(ps[i][t-1]+h*(qs[i+1][t-1]+qs[i-1][t-1]-2*qs[i][t-1]+b*((qs[i+1][t-1]-qs[i][t-1])**3-(qs[i][t-1]-qs[i-1][t-1])**3)))
        except OverflowError:
            print("** FATAL: OVERFLOW **")
            print("* diag: t =", t)
            # print(ps[i][t-1])
            # print(h*qs[i+1][t-1]+qs[i-1][t-1]-2*qs[i][t-1])
            # print(h*a*((qs[i+1][t-1]-qs[i][t-1])**2))
            # print(-h*a*((qs[i][t-1]-qs[i-1][t-1])**2))
            exit()
        #print("t=", t, "i=", i, qs[i+1][t-1], qs[i][t-1], h*a*((qs[i+1][t-1]-qs[i][t-1])**2))

        qs[i].append(qs[i][t-1] + h*(ps[i][t]))
    t = t + 1

# TRANSFORM BACK TO NORMAL MODES
Q = lambda k, t: (2/n)**(1/2) * sum([sin(pi*k*i/(n+1))*qs[i][t] for i in range(1, n+1)])
dQ = lambda k, t: (2/n)**(1/2) * sum([sin(pi*k*i/(n+1))*ps[i][t] for i in range(1, n+1)])

def E(k,t):
    return (0.5 * dQ(k,t)**2 + 0.5 * (omega(k)**2) * (Q(k,t)**2))

####################################################################
# PLOTTING ROUTINES
####################################################################
import matplotlib.pyplot as plot

t = 0
(E9, E10, E11, E12, E13) = ([], [], [], [], [])
Ts = [] # for pyplot plotting
sumE9 = 0
sumE10 = 0
sumE11 = 0
sumE12 = 0
sumE13 = 0
while t <= L:
    Ts.append(t*h/TU)
    E9.append(E(9,t))
    E10.append(E(10,t))
    E11.append(E(11,t))
    E12.append(E(12,t))
    E13.append(E(13,t))
    t = t + 1

# 2- problem information
# print(E1[159])
# print(E1[159], E10[159], E11[159], E12[159])

fig, ax = plot.subplots()
plot.plot(Ts, E9, label='E9', lw=1)
plot.plot(Ts, E10, label='E10', lw=1)
plot.plot(Ts, E11, label='E11', lw=1)
plot.plot(Ts, E12, label='E12', lw=1)
plot.plot(Ts, E13, label='E13', lw=1)
plot.legend()
ax.set_yscale('log')

plot.show()

# 5 - Q2=20
# for i in range(int(500*M/h), int(L)):
#     sumE1 = E1[i] + sumE1
#     sumE10 = E10[i] + sumE10
#     sumE11 = E11[i] + sumE11
#     sumE12 = E12[i] + sumE12

# avgE1 = sumE1/(-1000*TU/h+int(L))
# avgE10 = sumE10/(-1000*TU/h+int(L))
# avgE11 = sumE11/(-1000*TU/h+int(L))
# avgE12 = sumE12/(-1000*TU/h+int(L))

# print(avgE1, avgE10, avgE11, avgE12)