# Electric Field Configuration
LightMode = "elliptical"
if(LightMode == "linear"):
    Ex = lambda t: -0.075525 * (cos(0.00712509*t)**2) * cos(0.057*t) + 0.0188815 * cos(0.00712509*t) * sin(0.00712509*t) * sin(0.057*t)
    Ey = lambda t: 0
    E_max  = 0.0657197
elif(LightMode == "elliptical"):
    Ex = lambda t: -0.0675516 * (cos(0.00712509*t))**2 * cos(0.057*t) + 0.0168881 * sin(0.00712509*t) * sin(0.057*t) * cos(0.00712509*t)
    Ey = lambda t: 0.00844405 * cos(0.00712509*t) * cos(0.057*t) * sin(0.00712509*t) + 0.0337758 * (cos(0.00712509*t))**2 * sin(0.057*t)
    E_max = 0.0675516 # (t = 1.16205e-8)

# Potential
Ip =  0.5

# Initial/Final Time Configuration
T0 = -2 * 110.23
Tf =  2 * 110.23

# Sampling Rate
S  = 1000

# How many electron samples per time sample?
ES = 800

# Time-evolution timestep
TS = 1e-1

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

    return (x, y, vx, vy)

# [(4-tuple)] genElectronICMulti(t)
# Generate ES-many electron initial conditions for given time.
def genElectronICMulti(t):
    res = [(t, genElectronIC(t)) for i in range(ES)]
    res = [x for x in res if x[1] != None]

    return res

####################################################################
# TIME EVOLUTION OPERATORS
####################################################################
# (t, (x, y, vx, vy)) electronTimeEvolLaser((t, (x, y, vx, vy)))
# Evolves given electron through time 
def electronTimeEvolLaser(inpt):
    if(inpt[0] >= Tf):
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

    return (t, (x, y, vx, vy))

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

ICs = list(map(genElectronICMulti, Ts))

print(" OK", flush=True)

print("Flattening Data Structure...", end='', flush=True)
# Flatten list -- as they're 2,4-tuples with time built-in we no longer need to separate by
# initial condition
# Es = Electron Array
Es = sum(ICs, [])
print(" OK", flush=True)

print("Time Evolution (Laser T)...", end='', flush=True)

# Evolve all through time in multiprocessing
# Using timestep TS
Es = list(map(electronTimeEvolLaser, Es))

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

plot.hexbin(PInfXs, PInfYs, extent=(-1.5, 1.5, -1.5, 1.5), cmap="winter")
plot.colorbar()
plot.show()

# Generate Histograms for PInfXs & PInfYs
plot.hist(PInfXs, bins=150)
plot.show()

plot.hist(PInfYs, bins=150)
plot.show()
