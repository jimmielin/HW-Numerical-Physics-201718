####################################################################
# Computational Physics, 2017-18 Sem1
# HW-3 Ex-3-2b
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
# Now Playing: Top of the World - The Carpenters
####################################################################

dk = 0.01
ks = [-1.5 + i * dk for i in range(0, 301)]
kprojs = []

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
# For HW-3 Ex-3-2
#
# !!! MUST KEEP THIS THE SAME AS EX3-2
# THIS IS ONLY A UTILITY FILE FOR PLANE WAVE PROJECTION
####################################################################

# x range
xMin = -2000
xMax =  2000
dx   =     0.1

# t range
tMin = 0
tMax = 200 # 2Npi/w
dt   = 0.05  
NX   = int((xMax - xMin)/dx)
NT   = int((tMax - tMin)/dt)

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

# Plane wave e^(ikx)
for j in range(len(ks)):
	kvec = [e**(1j * ks[j] * (xMin + i * dx)) for i in range(NX)]
	kvec = vecNormalize(kvec)

	print(kvec[:5], kvec[(-5):])
	# kprojs.append(vecDotProduct(kvec, fvec))