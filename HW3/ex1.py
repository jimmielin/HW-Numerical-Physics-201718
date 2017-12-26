####################################################################
# Computational Physics, 2017-18 Sem1
# HW-3 Ex-1
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
# Now Playing: Love the Way You Lie - Eminem
####################################################################

####################################################################
# CONFIGURABLE PARAMETERS
# For HW-3 Ex-1
####################################################################

# Electric Field Configuration
LightMode = "linear"
if(LightMode == "linear"):
    Ex = lambda t: 0.075525 * (cos(0.00712509*t)**2) * cos(0.057*t) - 0.0188815 * cos(0.00712509*t) * sin(0.00712509*t) * sin(0.057*t)
    Ey = lambda t: 0
elif(LightMode == "elliptical"):
    Ex = lambda t: 0.0675516 * (cos(0.00712509*t))**2 * cos(0.057*t) - 0.0168881 * sin(0.00712509*t) * sin(0.057*t) * cos(0.00712509*t)
    Ey = lambda t: -0.00844405 * cos(0.00712509*t) * cos(0.057*t) * sin(0.00712509*t) - 0.0337758 * (cos(0.00712509*t))**2 * sin(0.057*t)
else:
    raise ValueError("LightMode must be either linear or elliptical")

# Potential
Ip =  0.5

# Initial/Final Time Configuration
T0 = -2 * 110.23
Tf =  2 * 110.23

# Sampling Rate
S  = 100

# How many electron samples per time sample?
ES = 100

# Time-evolution timestep
TS = 1e-3

# DIAGNOSTICS:
#   Set to 0 for quiet (production mode)
#   If > 0, various levels (100, 200, 1000) will return gradually
#   verbose diagnostic messages for the purpose of debugging.
# Not all functions support a diagnostic output.
debugLevel = 1000

# MULTITHREADING SUPPORT:
#   To speed up the program, enable multithreading below and configure
#   the appropriate core count (mtCoreCount)
#   - For debugging, turn off multithreading.
multiThreading = True
mtCoreCount    = 3
if(multiThreading):
    from multiprocessing import pool

####################################################################
# RANDOMIZER FUNCTIONS
# LCGRand() -> [0, 1] random number
####################################################################
# Provide a seed
LCGSeed = None

# If no seed, then seed with time
if(LCGSeed == None):
    import time
    LCGSeed = int(time.time())

# LCGRand()
# Generate random number in [0, 1] using the "Linear Congruence Generator"
# (1951, D.H. Lehmer)
#
# Using Parameters from Numerical Recipes:
# m = 2**32, a = 1664525, c = 1013904223
# But this is customizable.
#
# It is IMPERATIVE to seed LCGRand() with a truly random seed value,
# before starting its usage.
# For the sole purposes of HW-3, LCG is seeded with the system time,
# but you SHOULD NOT USE THIS FOR CRYPTOGRAPHIC PURPOSES
# (hell, you should not be using RNGs at all!)

# Internal State.
LCGInternalState = LCGSeed
def LCGRand():
    global LCGSeed, LCGInternalState

    # Numerical Recipes Parameters for LCG
    a = 1664525
    c = 1013904223
    m = 2**32

    LCGInternalState = (a * LCGInternalState + c) % m
    return LCGInternalState / m

# Do not remove.
# This runs the LCGRand a few times to remove the initial bias
# caused by time being a not-so-fast-moving variable.
for i in range(1, 255):
    LCGRand()

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
# DIFFERENTIAL EQUATION SOLVING FUNCTIONS (EULER)
# XY 2-DIMENSIONAL
# Mostly adapted from HW-1 Ex-5
####################################################################
# Euler Forward Method, for XY coupled situation
# eulerForwardSolveXY(\t, yx, yy -> f(t,yx), \t, yy, yx -> f(t,yy), a, b, x0, y0, dt = None)
# If timestep dt is none, then a default value will be used given N = 1000 partitions
# Returns the values for all the partitioned f(t_i) in a tuple (t_i, x_i, y_i)
#
# FIXME: default 1000 partitions is completely arbitrary
def eulerForwardSolveXY(fx, fy, a, b, x0, y0, dt = None):
    # Euler forward method
    k = 0
    if(dt != None):
        h = dt    # timestep given
    else:
        h = (b-a)/1000
    xn = a
    ys = [(a,x0,y0)]

    while((xn) <= b): # do a little nudging because python's double is a little too sensitive
        xn += h
        ys.append((xn, ys[-1][1] + h * fx(xn, ys[-1][1], ys[-1][2]), ys[-1][2] + h * fy(xn, ys[-1][2], ys[-1][1])))

    return ys


####################################################################
# INITIAL CONDITION GENERATOR
####################################################################
# (x, y, vx, vy) genElectronIC(t)
# Generate one random electron initial condition 4-tuple
# Positions x, y; velocities vx, vy
# Given a time variable.
# Uses acceptance-rejection sampling, which MAY BE SLOW...
def genElectronIC(t):
    _Ex = Ex(t)
    _Ey = Ey(t)
    E = (_Ex**2 + _Ey**2) ** (1/2)

    # Calculate angles to assign vx, vy later
    sinT = _Ey/E
    cosT = _Ex/E

    # Randomize a direction for v_perp
    vpD  = (1 if LCGRand() >= 0.5 else -1)

    # Distribution (without constant w.r.t. E-terms)
    Wd   = lambda vp: e**(-vp*vp/E)

    # Randomize v_perp's modulo (vpM) to fit distribution
    # WTarget = e**(-2/(3*E)) / (E*E) is even not necessary because
    # we only need to judge whether Target * \xi_2 <= f(\delta)
    # and this cancels out the regular terms. So we set target as 1.
    while(1 == 1):
        d = LCGRand() * 0.1    # Cap to a certain limit, noting v < c = 137
        ir_xi2 = LCGRand()
        ir_Wd = Wd(d)
        if(ir_xi2 <= ir_Wd):
            vp = d
            break

    # Now we have vp (v_perp) and its direction
    vp = vp * vpD

    # Compute x, ys. Direction is already implied, now use the fact that
    # (vx, vy) = (-vp*sinT, vp*cosT), the minus sign is exactly for the perpendicular
    vx = -vp*sinT
    vy =  vp*cosT

    # Compute position, which is *opposite to* the direction of E
    # and distance is dictated by |r0| = Ip/E (Ip = 0.5)
    r0 = Ip/E
    x  = -r0*cosT
    y  = -r0*sinT

    return (x, y, vx, vy)

# [(4-tuple)] genElectronICMulti(t)
# Generate ES-many electron initial conditions for given time.
def genElectronICMulti(t):
    if(debugLevel >= 500):
        print("* diag genElectronICMulti: generating samples for t =", t, flush=True)
    return [(t, genElectronIC(t)) for i in range(ES)]

####################################################################
# TIME EVOLUTION OPERATORS
####################################################################
# (t, (x, y, vx, vy)) electronTimeEvolLaser((t, (x, y, vx, vy)))
# Evolves given electron through time 
def electronTimeEvolLaser(inpt):
    if(debugLevel >= 1000):
        print("* diag electronTimeEvolLaser: input ", inpt, flush=True)

    if(inpt[0] >= Tf):
        # already fully evolved
        return inpt

    # To evolve through time, first integrate acceleration (dv/dt) -> velocity...
    # Force generators:
    # EM field --
    _Ex = lambda t, x, y: Ex(t)
    _Ey = lambda t, x, y: Ey(t)
    # Coulomb Potential (in force form 1/r^2) -- 
    _Cmod = lambda t, x, y: 1/(x**2 + y**2 + 0.04)

    # Initialize Internal State
    try:
        x = inpt[1][0]
        y = inpt[1][1]
        vx = inpt[1][2]
        vy = inpt[1][3]
    except IndexError:
        print(inpt)
        raise IndexError("** Above are diagnostics for IndexError **")

    # Number of samples
    NS = int(((Tf - inpt[0])/TS)//1 + 1)
    Ts = [inpt[0] + (Tf - inpt[0])/NS * i for i in range(NS + 1)]

    # Evolve for each timestep:
    for t_ptr in range(len(Ts)):
        t = Ts[t_ptr]
        r = (x**2 + y**2) ** (1/2)

        # Generate forces...
        Cmod = _Cmod(t, x, y)
        Fx = _Ex(t, x, y) - x/r * Cmod
        Fy = _Ey(t, x, y) - y/r * Cmod

        # Integrate one step towards F-direction to get velocities...
        (vx, vy) = (vx + Fx * TS, vy + Fy * TS)

        # Integrate one step towards V-direction to get Xs
        (x, y) = (x + vx * TS, y + vy * TS)

        # ...debug?
        if(debugLevel >= 1000):
            print("* diag electronTimeEvolLaser evolved under timestep for ", (t, x, y, vx, vy), flush=True)

    return (t, (x, y, vx, vy))


# Main Runner
# This clause is ESSENTIAL for Multiprocessing to perform properly.
# Otherwise threads will not fork correctly.
if __name__ == '__main__':
    # GENERATE INITIAL SAMPLES
    # Perform a sampling of T
    Ts = [T0 + (Tf - T0)/S * i for i in range(S)]

    # Create a list to store samples of electron INITIAL CONDITIONS in memory,
    # in the following format: keys are ordered in the same spacing as Ts
    # e.g. Ts[0] = -220.46s
    # and ICs[0] = [ ... (t, (x, y, vx, vy)) ... ] in a list of 2,4-tuples
    # length is same as len(Ts)
    ICs = [[]] * len(Ts)

    if(multiThreading == True):
        MPool = pool.Pool(mtCoreCount)
        ICs   = MPool.map(genElectronICMulti, Ts)
    else:
        raise NotImplementedError("Non-multiprocessing based code is unsupported at the moment.")

    # Flatten list -- as they're 2,4-tuples with time built-in we no longer need to separate by
    # initial condition
    # Es = Electron Array
    Es = sum(ICs, [])

    # Evolve all through time in multiprocessing
    # Using timestep TS
    if(multiThreading == True):
        # We provide a list of 2,4-tuples where (t, (x, y, vx, vy))
        # this is so we keep track of each electron in time and facilitates parrying
        # in multiprocessing.pool.
        MPool = pool.Pool(mtCoreCount)
        Es = MPool.map(electronTimeEvolLaser, Es)

        print(Es)
    else:
        raise NotImplementedError("Non-multiprocessing based code is unsupported at the moment.")


