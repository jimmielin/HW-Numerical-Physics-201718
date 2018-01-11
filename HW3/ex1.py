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
LightMode = "elliptical"
if(LightMode == "linear"):
    Ex = lambda t: -0.075525 * (cos(0.00712509*t)**2) * cos(0.057*t) + 0.0188815 * cos(0.00712509*t) * sin(0.00712509*t) * sin(0.057*t)
    Ey = lambda t: 0
    # Ex_max = 0.0419221 (t = +/- 102.851)
    E_max  = 0.0657197
elif(LightMode == "elliptical"):
    Ex = lambda t: -0.0675516 * (cos(0.00712509*t))**2 * cos(0.057*t) + 0.0168881 * sin(0.00712509*t) * sin(0.057*t) * cos(0.00712509*t)
    Ey = lambda t: 0.00844405 * cos(0.00712509*t) * cos(0.057*t) * sin(0.00712509*t) + 0.0337758 * (cos(0.00712509*t))**2 * sin(0.057*t)
    # Ex_max = 0.0374962 (t = +/- 102.851)
    # Ey_max = 0.0326369 (t = +/- 25.9684)
    E_max = 0.0675516 # (t = 1.16205e-8)
else:
    raise ValueError("LightMode must be either linear or elliptical")

# Potential
Ip =  0.5

# Initial/Final Time Configuration
T0 = -2 * 110.23
Tf =  2 * 110.23

# Sampling Rate
S  = 8000

# How many electron samples per time sample?
ES = 4000

# Time-evolution timestep
TS = 1e-1

# DIAGNOSTICS:
#   Set to 0 for quiet (production mode)
#   If > 0, various levels (100, 200, 1000) will return gradually
#   verbose diagnostic messages for the purpose of debugging.
# Not all functions support a diagnostic output.
debugLevel = 0

# MULTITHREADING SUPPORT:
#   To speed up the program, enable multithreading below and configure
#   the appropriate core count (mtCoreCount)
#   - For debugging, turn off multithreading.
multiThreading = True
mtCoreCount    = 48
if(multiThreading):
    from multiprocessing import pool

# GENERATE PLOTS?
#   If yes, matplotlib.pyplot must be provided.
genPlots      = True
if(genPlots):
    import matplotlib.pyplot as plot

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

    # Distribution (WITH constant w.r.t. E-terms because it affects the overall E-vprep distribution)
    Wd   = lambda vp: 1/(E**2) * e**(-2/(3*E)) * e**(-vp*vp/E)
    WdM  = 1/(E_max**2) * e**(-2/(3*E_max))# * e**(-vp*vp/E_max)

    # Randomize v_perp's modulo (vpM) to fit distribution
    # WTarget = e**(-2/(3*E)) / (E*E) is even not necessary because
    # we only need to judge whether Target * \xi_2 <= f(\delta)
    # and this cancels out the regular terms. So we set target as 1.
    while(1 == 1):
        ir_xi1 = LCGRand() * 10    # Cap to a certain limit, noting v < c = 137
        ir_xi2 = WdM * LCGRand()
        ir_Wd = Wd(ir_xi1)
        if(ir_xi2 <= ir_Wd):
            vp = ir_xi1
            break
        else:
            return None

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

    # print(r0)

    # return (x, y, vx, vy, vp*vp/E) # for debug only
    return (x, y, vx, vy)

# [(4-tuple)] genElectronICMulti(t)
# Generate ES-many electron initial conditions for given time.
def genElectronICMulti(t):
    if(debugLevel >= 500):
        print("* diag genElectronICMulti: generating samples for t =", t, flush=True)
    res = [(t, genElectronIC(t)) for i in range(ES)]
    # If samples are unavailable, remove these elements...
    res = [x for x in res if x[1] != None]
    
    if(debugLevel >= 100):
       print("* diag genElectronICMulti: generated samples for t, num =", t, len(res), flush=True)

    return res

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
    x = inpt[1][0]
    y = inpt[1][1]
    vx = inpt[1][2]
    vy = inpt[1][3]

    # Number of time samples
    NS = int(((Tf - inpt[0])/TS)//1 + 1)
    Ts = [inpt[0] + (Tf - inpt[0])/NS * i for i in range(NS + 1)]

    # Evolve for each timestep:
    for t_ptr in range(len(Ts)):
        t = Ts[t_ptr]
        r = (x**2 + y**2) ** (1/2)

        # Generate forces...
        Cmod = _Cmod(t, x, y)
        Fx = - _Ex(t, x, y) - x/r * Cmod
        Fy = - _Ey(t, x, y) - y/r * Cmod

        # Integrate one step towards F-direction to get velocities...
        vx = vx + Fx * TS
        vy = vy + Fy * TS

        # Integrate one step towards V-direction to get Xs
        x = x + vx * TS 
        y = y + vy * TS

        # (x, y, vx, vy) = (x + vx * TS, y + vy * TS, vx + Fx * TS, vy + Fy * TS)

        # ...debug?
        if(debugLevel >= 1000):
            print("* diag electronTimeEvolLaser evolved under timestep for ", (t, x, y, vx, vy, Fx, Fy), flush=True)

    if(debugLevel > 0):
        print("*", flush=True, end='')
    return (t, (x, y, vx, vy))


# Main Runner
# This clause is ESSENTIAL for Multiprocessing to perform properly.
# Otherwise threads will not fork correctly.
if __name__ == '__main__':
    print("********************************************")
    print("* COMPUTATIONAL PHYSICS, HW-3 Ex-1         *")
    print("* (c) 2017-2018 Haipeng Lin                *")
    print("********************************************")
    if(multiThreading == True):
        print("* MULTITHREADING IS ENABLED.               *")
    else:
        print("* [WARNING] Running on SINGLE-THREADED MODE*")

    print("********************************************\n")

    print("Generating initial samples...", end='', flush=True)
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

    print(" OK", flush=True)

    # DEBUG -- Statistical Counting in relation to E
    # ICCounts = [len(ICss) for ICss in ICs]
    # import matplotlib.pyplot as plt
    # fig, ax = plt.subplots()
    # Es = [(Ex(t)**2 + Ey(t)**2)**(1/2) for t in Ts]
    # ax.scatter(Es, ICCounts, s=1)
    # plt.show()

    print("Flattening Data Structure...", end='', flush=True)
    # Flatten list -- as they're 2,4-tuples with time built-in we no longer need to separate by
    # initial condition
    # Es = Electron Array
    Es = sum(ICs, [])
    print(" OK", flush=True)

    # DEBUG -- Statistical Counting in relation to vp/E (only when you have a 5th-dimension argument in 2,5-tuple)
    # vp2dE = [E[1][4] for E in Es]
    # import matplotlib.pyplot as plt
    # fig, ax = plt.subplots()
    # ax.hist(vp2dE, bins = 100)
    # plt.show()

    # DEBUG -- View final positions of electrons (flattened structure)
    # xs = []
    # ys = []
    # vxs = []
    # vys = []

    # for i in range(len(Es)):
    #     xs.append(Es[i][1][0])
    #     ys.append(Es[i][1][1])
    #     vxs.append(Es[i][1][2])
    #     vys.append(Es[i][1][3])

    # plot.hexbin(xs, ys, cmap="hot")
    # plot.colorbar()
    # plot.show()

    # plot.hexbin(vxs, vys, extent=(-1.5, 1.5, -1.5, 1.5), cmap="hot")
    # plot.colorbar()
    # plot.show()

    # exit()

    print("Time Evolution (Laser T)...", end='', flush=True)

    # Evolve all through time in multiprocessing
    # Using timestep TS
    if(multiThreading == True):
        # We provide a list of 2,4-tuples where (t, (x, y, vx, vy))
        # this is so we keep track of each electron in time and facilitates parrying
        # in multiprocessing.pool.
        MPool = pool.Pool(mtCoreCount)
        Es = MPool.map(electronTimeEvolLaser, Es)
    else:
        raise NotImplementedError("Non-multiprocessing based code is unsupported at the moment.")

    print(" OK", flush=True)
    print("Flattening Structure for p_\\infty...", end='', flush=True)
    
    #   0       10 11  12  13
    # [(220.46, (x, y, vx, vy))]
    # Non-Laser Time Evolution can be algebraically solved for (see HW-3 Ex-1 eqs. (2))
    # Iterate through electron samples (Es) and resolve for their PInfX/PInfYs
    PInfXs = []
    PInfYs = []

    for i in range(len(Es)):
        x = Es[i][1][0]
        y = Es[i][1][1]
        vx = Es[i][1][2]
        vy = Es[i][1][3]

        pfSq = vx ** 2 + vy ** 2
        rf   = (x ** 2 + y ** 2) ** (1/2)
        Lz   = (x * vy - y * vx)
        LzSq = Lz ** 2

        PInfSq = pfSq - 2/rf
        PInfM = PInfSq ** (1/2)

        # If total energy < 0, then it hasn't left the atom but instead has been excited to another state
        if(PInfSq < 0):
            continue

        PInfX = PInfM / (1 + PInfSq * LzSq) * (PInfM * (LzSq * vx + y * Lz / rf) - Lz * vy + x / rf)
        PInfY = -PInfM / (1 + PInfSq * LzSq) * (PInfM * (LzSq * vy - x * Lz / rf) + Lz * vx + y / rf)

        # print(vx, PInfX, vy, PInfY) # not too much difference either...

        PInfXs.append(PInfX)
        PInfYs.append(PInfY)

    print(" OK", flush=True)

    if(genPlots == True):
        plot.hexbin(PInfXs, PInfYs, extent=(-1.5, 1.5, -1.5, 1.5), cmap="plasma")
        plot.colorbar()
        plot.show()

        # Generate Histograms for PInfXs & PInfYs
        histX = plot.hist(PInfXs, bins=150)
        plot.show()

        plot.hist(PInfYs, bins=150)
        plot.show()
    else:
        print(EPInfsPlot)
