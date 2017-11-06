# $ RELEASE $
# $ 201711060205Z $ rel01
# $ Signed-Off-By: Haipeng Lin <jimmie.lin@gmail.com>
####################################################################
# Computational Physics, 2017-18 Sem1
# HW-1 Ex-7
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
# This program is Python 3 (3.6.2 on darwin)
# compatible, and does not support Python 2.
#
# Now Playing: Forgettable - Project 46
#              喜欢恋爱 - 卢巧音
####################################################################

n = 10000

# Initialize the matrices *using MATHEMATICAL INDEXES*
# This is to prevent index confusion.

# Primary Diagonal & X Solution Array
M = [6] * (n + 1)
M[0] = 0 # padding index
M[1] = 5
M[n] = 5

# Secondary Diagonal
A = [4] * (n + 1)
A[0] = 0 # non-significant padding index
A[1] = 0 # non-significant padding index

# Tertiary Diagonal
B = [1] * (n + 1)
B[0] = 0 # non-significant padding index
B[1] = 0 # non-significant padding index
B[2] = 0 # non-significant padding index

# Perform first step calculation
M[1] = 60 / M[1]
for i in range(2, n):
    # B[i]X[i-2]+A[i]X[i-1]+M[i]X[i]=120
    # X[i]=1/M[i] * (120-B[i]X[i-2]-A[i]X[i-1])
    M[i] = 1/M[i] * (120 - B[i] * M[i-2] - A[i] * M[i-1])
M[n] = 1/M[n] * (60 - B[n] * M[n-2] - A[n] * M[n-1])

print(M[1:])