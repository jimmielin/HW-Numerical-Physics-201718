####################################################################
# Computational Physics, 2017-18 Sem1
# HW-2 Ex-2(3)
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
# This program is Python 3 (3.6.2 on MSYS_NT)
# compatible, and does not support Python 2.
#
# Now Playing: Kids In The Dark - All Time Low
####################################################################

from pylitematrix.pylitematrix import Matrix

# Configuration Parameters
# Use Hermite N: 96, 48, 24, 12 are supported ONLY
N = 12

# Configure constant factor s (in ex.2-3 values are 0.6, 1.0)
s = 1

# Configure potential V(z)
V = lambda z: z**2 + z**4

## !!! FOR AN EXPLANATION OF WHY THESE SOLUTIONS ARE BUILT IN AND HOW THEY WERE GENERATED, SEE DOCUMENTATION !!! ##
## !!! THE ROOTS ARE PRE-SOLVED USING MY PROGRAM !!! ##
# Built-in Solutions for Hermite Polynomials
# So we don't spend time solving them all over again...
H96 = [-13.116130021662878, -12.529230746683542, -12.044806741134922, -11.613620828842913, -11.216875194202162, -10.844888073853587, -10.491854539276538, -10.153951274713688, -9.828492746191714, -9.513501769412425, -9.20746935337176, -8.909210581459599, -8.617773244214746, -8.33237730717006, -8.052373330719712, -7.777213031061239, -7.506427894427756, -7.239613294008353, -6.976416464285181, -6.716527240496519, -6.459670819555122, -6.20560202471831, -5.954100706411726, -5.704968013525535, -5.4580233400707865, -5.2131018018189605, -4.970052133166773, -4.728734920352931, -4.489021106216995, -4.250790715903187, -4.0139317636282685, -3.7783393087968626, -3.5439146360273117, -3.310564538524305, -3.0782006880468185, -2.846739077724783, -2.616099526362389, -2.38620523476978, -2.1569823861931283, -1.9283597841455025, -1.7002685219383507, -1.4726416790227095, -1.2454140399065081, -1.0185218319476799, -0.7919024787540168, -0.5654943662667489, -0.33923661887888523, -0.11306888315145194, 0.11306888315145194, 0.33923661887888523, 0.5654943662667489, 0.7919024787540168, 1.0185218319476799, 1.2454140399065081, 1.4726416790227095, 1.7002685219383507, 1.9283597841455025, 2.1569823861931283, 2.38620523476978, 2.616099526362389, 2.846739077724783, 3.0782006880468185, 3.310564538524305, 3.5439146360273117, 3.7783393087968626, 4.0139317636282685, 4.250790715903187, 4.489021106216995, 4.728734920352931, 4.970052133166773, 5.2131018018189605, 5.4580233400707865, 5.704968013525535, 5.954100706411726, 6.20560202471831, 6.459670819555122, 6.716527240496519, 6.976416464285181, 7.239613294008353, 7.506427894427756, 7.777213031061239, 8.052373330719712, 8.33237730717006, 8.617773244214746, 8.909210581459599, 9.20746935337176, 9.513501769412425, 9.828492746191714, 10.153951274713688, 10.491854539276538, 10.844888073853587, 11.216875194202162, 11.613620828842913, 12.044806741134922, 12.529230746683542, 13.116130021662878]

H48 = [-8.975315081931686, -8.310752190704784, -7.759295519765774, -7.26604655416435, -6.810064578074141, -6.3805640961864105, -5.971072225013545, -5.577316981223729, -5.196287718792364, -4.8257572281332095, -4.464014546934459, -4.109704603560591, -3.761726490228358, -3.4191659693638847, -3.0812489886451058, -2.747308624822383, -2.4167609048732164, -2.0890866609442766, -1.7638175798953, -1.4405252201375651, -1.1188121524021566, -0.7983046277785623, -0.4786463375944961, -0.15949293584886248, 0.15949293584886248, 0.4786463375944961, 0.7983046277785623, 1.1188121524021566, 1.4405252201375651, 1.7638175798953, 2.0890866609442766, 2.4167609048732164, 2.747308624822383, 3.0812489886451058, 3.4191659693638847, 3.761726490228358, 4.109704603560591, 4.464014546934459, 4.8257572281332095, 5.196287718792364, 5.577316981223729, 5.971072225013545, 6.3805640961864105, 6.810064578074141, 7.26604655416435, 7.759295519765774, 8.310752190704784, 8.975315081931686]

H24 = [-6.01592556142574, -5.259382927668044, -4.625662756423787, -4.05366440244815, -3.5200068130345246, -3.0125461375655647, -2.5238810170114268, -2.049003573661699, -1.5842500109616942, -1.1267608176112451, -0.6741711070372123, -0.22441454747251557, 0.22441454747251557, 0.6741711070372123, 1.1267608176112451, 1.5842500109616942, 2.049003573661699, 2.5238810170114268, 3.0125461375655647, 3.5200068130345246, 4.05366440244815, 4.625662756423787, 5.259382927668044, 6.01592556142574]

H12 = [-3.889724897869782, -3.0206370251208896, -2.2795070805010598, -1.5976826351526048, -0.9477883912401638, -0.31424037625435913, 0.31424037625435913, 0.9477883912401638, 1.5976826351526048, 2.2795070805010598, 3.0206370251208896, 3.889724897869782]

# Choose the Hs solution list for what we want to do
# Also padd the list so there are mathematical indeces instead
if(N == 96):
	Hs = [None] + H96
elif(N == 48):
	Hs = [None] + H48
elif(N == 24):
	Hs = [None] + H24
elif(N == 12):
	Hs = [None] + H12
else:
	raise ValueError("Unsupported N -- provide HXX?")

mp = Matrix(1, 1, lambda i, j: 1)

# Generate Matrices for solving the Eigenvalue Problem.
matrixT = Matrix(N, N, lambda i, j: ((4*N - 1 - 2*Hs[i]*Hs[i])/6 if i == j else ((1 if (i + j) % 2 == 0 else -1)*(2/(Hs[i]-Hs[j])**2 - 0.5)))/(s*s))

matrixV = Matrix(N, N, lambda i, j: V(Hs[i]*s) if i == j else 0)

# For a vector, provide the norm (2-norm)
# If you pass a not 1x (x1) matrix, the behavior is undefined
def norm(v):
	return sum([x**2 for x in (v.toList() if isinstance(v, Matrix) else v)]) ** (1/2)

# For a vector, provide the norm (infty-norm)
# If you pass a not 1x (x1) matrix, the behavior is undefined
def normInfty(v):
	return max((v.toList() if isinstance(v, Matrix) else v))

# Householder Transformation
# householder(v) = H
# Per Zhou et, al. Tsinghua University Press.
#
# If you input v as Matrix, it will be checked for vector
# and transformed as a list. But you can also just provide the list.
def householder(v):
	if(not isinstance(v, list)):
		if(v.nrows() != 1 and v.ncols() != 1):
			raise ValueError("Input for householder must be a vector!")

		# What dimensions are we? Use lists from now on, Matrixes are tedious
		# and slow when operating on 1x x1 vectors
		V = v.toList()
	else:
		V = v

	n = len(V)

	# Normalize vector v
	vM = normInfty(v)
	if(abs(vM) < 1e-10):
		# Too small -- we can consider this zero and return a 1 transformation
		return mp.id_(n)

	nV = [i/vM for i in V]

	# Padd so there are mathematical indexes...
	nV = [None] + nV

	sigma = sum([i**2 for i in nV[2:]])
	hV = [0, 1]
	hV = hV + nV[2:]
	if(sigma == 0):
		beta = 0
	else:
		try:
			alpha = (nV[1]**2 + sigma)**(1/2)
		except OverflowError:
			print("\n\n***************************************")
			print("A OVERFLOWERROR OCCURRED IN THE HOUSEHOLDER TRANSFORMATION ROUTINE.")
			print("Debugging information is available below")
			print("***************************************")
			print("vM (infty norm) =", vM)
			print("nV[1] =", nV[1])
			print("sigma =", sigma)
			raise ValueError("Result too large... aborting with error tracelog")
		if(nV[1] <= 0):
			hV[1] = nV[1] - alpha
		else:
			hV[1] = -sigma/(nV[1] + alpha)
		beta = 2*(hV[1]**2)/(sigma + hV[1]**2)
	hV = [x/hV[1] for x in hV]
	hV = hV[1:]

	hVm = mp.fromList_(n, 1, hV)
	hVmT = hVm.transpose()

	bMapd = (hVm * hVmT).map(lambda i: -i * beta)

	H = Matrix(n, n, lambda i, j: 1 if i == j else 0) + bMapd

	return H

# (Matrix Q, Matrix R) qrDecomposition(Matrix A)
# Performs QR Decomposition of Matrix A, iteratively over every submatrix dimension
def qrDecomposition(A):
	# print("A:\n", A, "\n\n")
	cols = A.ncols()
	rows = A.nrows()
	HHs = []
	HHsMult = mp.id_(rows)
	for i in range(1, cols + 1):
		# Get the sub-rowvec using subMatrix(sr, er, sc, ec)
		subvec = A.subMatrix(i, rows, i, i)
		# print(subvec)
		# Perform Householder Transformation
		hh = householder(subvec)

		# Padd hh to the right dimension (should be rows x rows)
		tl = mp.id_(i - 1)
		tr = mp.zeros_(i - 1, rows - i + 1)

		# print("tl:\n", tl)
		# print("tr:\n", tr)

		bl = mp.zeros_(rows - i + 1, i - 1)
		br = hh

		# print("bl:\n", bl)
		# print("br:\n", br)

		pHH = mp.joinBlocks_(tl, tr, bl, br)

		# print("pHH:\n", pHH)

		HHsMult = HHsMult * pHH

		# Modify A to be a different A...
		A = pHH * A
		# print("A':\n", A, "\n\n")
		
	Q = HHsMult
	R = A
	return (Q, R)

# (Matrix H) hessenbergMatrix(Matrix A)
# Decompose Matrix A into a upper-Hessenberg Matrix for an improved QR method using
# Householder/Givens transformations
def hessenbergMatrix(A):
	cols = A.ncols()
	rows = A.nrows()
	for i in range(1, cols): # NOT cols + 1, see Tsinghua, p. 156
		# Get the sub-rowvec using subMatrix(sr, er, sc, ec)
		# note we only transform this to tridiag, so we start 1 row after
		subvec = A.subMatrix(i + 1, rows, i, i)

		# Perform Householder Transformation
		hh = householder(subvec)

		# Padd hh to the right dimension (should be rows x rows)
		tl = mp.id_(i)
		tr = mp.zeros_(i, rows - i)

		# print("tl:\n", tl)
		# print("tr:\n", tr)

		bl = mp.zeros_(rows - i, i)
		br = hh

		# print("bl:\n", bl)
		# print("br:\n", br)

		pHH = mp.joinBlocks_(tl, tr, bl, br)

		# Modify A to be a different A...
		A = pHH * A * pHH


	# # Perform a pass-filter to remove < e-10
	# for i in range(len(A.internal)):
	# 	for j in range(len(A.internal[i])):
	# 		if(A.internal[i][j] < 1e-10):
	# 			A.internal[i][j] = 0

	return A

# (List eigenvalues) qrSolveEigenvalues(Matrix A, Int MaxIterations)
# Solve Eigenvalues (returning a almost diagonal matrix) using QR Iteration Method
def qrSolveEigenvalues(A, MaxIterations = 100, Eigenvectors = False):
	iterations = 0
	L = A
	if(Eigenvectors):
		eVs = mp.id_(n)

	while(iterations < MaxIterations):
		(Q, R) = qrDecomposition(L)
		L = R*Q
		iterations += 1

		if(Eigenvectors):
			eVs = Q * eVs

		# print("*", end='', flush=True)
		print([L.internal[i][i] for i in range(0, L.nrows())], "\n")

	print('\n')

	# Get diagonal elements...
	if(Eigenvectors):
		return ([L.internal[i][i] for i in range(0, n)], eVs)
	else:
		return [L.internal[i][i] for i in range(0, n)]

# (List eigenvalues) qrSolveEigenvaluesWithTransform(Matrix A, Int MaxIterations)
# Solve Eigenvalues (returning a almost diagonal matrix) using QR Iteration Method,
# with preparation of a upper Hessenberg matrix for improved computational speed.
def qrSolveEigenvaluesWithTransform(A, MaxIterations = 100):
	iterations = 0
	L = hessenbergMatrix(A)
	while(iterations < MaxIterations):
		(Q, R) = qrDecomposition(L)
		L = R*Q
		iterations += 1


		# print("*", end='', flush=True)
		print([L.internal[i][i] for i in range(0, L.nrows())], "\n")

	print('\n')

	# Get diagonal elements...
	return [L.internal[i][i] for i in range(0, L.nrows())]


# (List eigenvalues) qrOffsetSolveEigenvalues(Matrix A, Int MaxIterations)
# Solve Eigenvalues (returning a almost diagonal matrix) using QR Iteration Method,
# with improved speed by introducing a offset to the iterations.
#
# Provide a Offset that should be preferably similar to an eigenvalue (known)
# of A. If an Offset is not provided, we'll use the bottom-rightmost "eigenvalue"
#
# How many iterations w/ offset?
# If Eigenvectors = True, then the eigenvectors are also returned, a tuple (Eigenvalues, Matrix Eigenvectors). There is a slight perf. hit.
def qrOffsetSolveEigenvalues(A, Offset = None, MaxIterations = 100, OffsetIterations = None, Eigenvectors = False):
	iterations = 0
	L = A
	n = A.nrows()

	if(Eigenvectors):
		eVs = mp.id_(n)

	while(iterations < MaxIterations):
		if(OffsetIterations == None or iterations < OffsetIterations):
			if(Offset != None):
				oM = mp.diag_(n, [Offset] * n)
				oMI = mp.diag_(n, [-Offset] * n)

			(Q, R) = qrDecomposition((L if Offset == None else L + oMI))
			L = (R*Q if Offset == None else R*Q + oM)
			irs = [L.internal[i][i] for i in range(0, n)] # intermediate result

			if(Offset == None):
				Offset = irs[-1] # use the last element
		else:
			if(Offset != None):
				oM = mp.diag_(n, [Offset] * n)
				oMI = mp.diag_(n, [-Offset] * n)

			(Q, R) = qrDecomposition((L if Offset == None else L + oMI))
			L = (R*Q if Offset == None else R*Q + oM)
			irs = [L.internal[i][i] for i in range(0, n)] # intermediate result

			if(Offset == None):
				Offset = irs[-1] # use the last element

		if(Eigenvectors):
			eVs = Q * eVs

		# print(irs, "\n")
		iterations += 1

	# print('\n')

	# Get diagonal elements...
	if(Eigenvectors):
		return ([L.internal[i][i] for i in range(0, n)], eVs)
	else:
		return [L.internal[i][i] for i in range(0, n)]

def qrOffsetSolveEigenvaluesWithTransform(A, Offset = None, MaxIterations = 100, OffsetIterations = None, Eigenvectors = False):
	iterations = 0
	L = hessenbergMatrix(A)
	n = A.nrows()

	if(Eigenvectors):
		eVs = mp.id_(n)

	while(iterations < MaxIterations):
		if(OffsetIterations == None or iterations < OffsetIterations):
			if(Offset != None):
				oM = mp.diag_(n, [Offset] * n)
				oMI = mp.diag_(n, [-Offset] * n)

			(Q, R) = qrDecomposition((L if Offset == None else L + oMI))
			L = (R*Q if Offset == None else R*Q + oM)
			irs = [L.internal[i][i] for i in range(0, n)] # intermediate result

			if(Offset == None):
				Offset = irs[-1] # use the last element
		else:
			if(Offset != None):
				oM = mp.diag_(n, [Offset] * n)
				oMI = mp.diag_(n, [-Offset] * n)

			(Q, R) = qrDecomposition((L if Offset == None else L + oMI))
			L = (R*Q if Offset == None else R*Q + oM)
			irs = [L.internal[i][i] for i in range(0, n)] # intermediate result

			if(Offset == None):
				Offset = irs[-1] # use the last element

		if(Eigenvectors):
			eVs = Q * eVs

		# print(irs, "\n")
		iterations += 1

	# print('\n')

	# Get diagonal elements...
	if(Eigenvectors):
		return ([L.internal[i][i] for i in range(0, n)], eVs)
	else:
		return [L.internal[i][i] for i in range(0, n)]

# Now we have to solve the Eigenvalue Problem
# (T+V) * Psi = E * Psi
matrixTpV = matrixT + matrixV
# print(matrixTpV.internal)
# print(qrSolveEigenvalues(matrixTpV))
# eigenvalues = (qrOffsetSolveEigenvalues(matrixTpV, Offset = 4.648812, MaxIterations = 100))
eigenvalues = qrSolveEigenvaluesWithTransform(matrixTpV, MaxIterations = 100)
eigenvalues.sort()
for i in range(len(eigenvalues)):
	print(i+1, ": ", eigenvalues[i])



###############  Testing code below, can safely ignore, use for debug ##################

# Test Householder Transformations
#print(householder(mp.fromList_(4, 1, [4, 3, 2, 1])))
#print(householder([4,3,2,1])*mp.fromList_(4,1,[4,3,2,1]))

# Test QR Decomposition
# A = mp.fromList_(5,5,[i for i in H96])
# print("A:\n", A)
# qrDecomp = qrDecomposition(A)
# print("Q:\n", qrDecomp[0]) # Q
# print('\n')
# print("R:\n", qrDecomp[1]) # R
# print('\n')
# print(qrDecomp[0] * qrDecomp[1]) # QR = A

# Test matrix prototype joinBlocks
# a = mp.fromList_(2, 2, [1,2,3,4])
# b = mp.fromList_(2, 3, [9,8,7,6,5,4])
# c = mp.fromList_(1, 2, [11, 13])
# d = mp.fromList_(1, 3, [89,64,0])

# Test Eigenvalues using QR Solving (Naive)
# Test data from Tsinghua, p. 154
# print(qrSolveEigenvalues(mp.fromList_(4, 4, [2.9766,0.3945,0.4198,1.1159,0.3945,2.7328,-0.3097,0.1129,0.4198,-0.3097,2.5675,0.6079,1.1159,0.1129,0.6079,1.7231])))


# print(mp.joinBlocks_(a,b,c,d))