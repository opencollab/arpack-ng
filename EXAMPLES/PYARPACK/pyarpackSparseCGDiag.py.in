#!/usr/bin/env python

import numpy as np
from pyarpack import sparseCGDiag as pyarpackSlv

# Build laplacian.

n = 4
i = np.array([], dtype='@PYINT@')
j = np.array([], dtype='@PYINT@')
Aij = np.array([], dtype='float64')
for k in range(n):
    for l in [k-1, k, k+1]:
        if l < 0 or l > n-1:
            continue
        i = np.append(i, np.@PYINT@(k)) # Casting value on append is MANDATORY or C++ won't get the expected type.
        j = np.append(j, np.@PYINT@(l)) # Casting value on append is MANDATORY or C++ won't get the expected type.
        if l == k:
            Aij = np.append(Aij, np.float64( 200.)) # Casting value on append is MANDATORY or C++ won't get the expected type.
        if l == k-1 or l == k+1:
            Aij = np.append(Aij, np.float64(-100.)) # Casting value on append is MANDATORY or C++ won't get the expected type.
for k, l, Akl in zip(i, j, Aij):
    print("A[", k, ",", l, "] =", Akl)
A = (n, i, j, Aij) # coo format: dimension, i 0-based indices, j 0-based indices, Aij values.

# Get and tune arpack solver.

arpackSlv = pyarpackSlv.double() # Caution: double <=> np.array(..., dtype='float64')
arpackSlv.verbose = 3 # Set to 0 to get a quiet solve.
arpackSlv.debug = 1 # Set to 0 to get a quiet solve.
arpackSlv.nbEV = 1
arpackSlv.nbCV = 2*arpackSlv.nbEV + 1
arpackSlv.mag = 'LM'
arpackSlv.maxIt = 200
arpackSlv.slvTol = 1.e-6
arpackSlv.slvMaxIt = 100

# Solve eigen problem.

rc = arpackSlv.solve(A)
assert rc == 0, "bad solve"
rc = arpackSlv.checkEigVec(A)
assert rc == 0, "bad checkEigVec"

# Print out results (mode selected, eigen vectors, eigen values, ...).

assert arpackSlv.nbEV == len(arpackSlv.val), "bad result"
print("\nresults:\n")
print("mode selected:", arpackSlv.mode)
print("nb iterations:", arpackSlv.nbIt)
print("Reverse Communication Interface time:", arpackSlv.rciTime, "s")
for val, vec in zip(arpackSlv.val, arpackSlv.vec):
    print("eigen value:", val)
    print("eigen vector:")
    print(vec)

#######################################################################################
print("\n##########################################################################\n")
#######################################################################################

# Build laplacian.

n = 8
i = np.array([], dtype='@PYINT@')
j = np.array([], dtype='@PYINT@')
Aij = np.array([], dtype='float32')
Bij = np.array([], dtype='float32')
for k in range(n):
    for l in [k-1, k, k+1]:
        if l < 0 or l > n-1:
            continue
        i = np.append(i, np.@PYINT@(k+1)) # Casting value on append is MANDATORY or C++ won't get the expected type.
        j = np.append(j, np.@PYINT@(l+1)) # Casting value on append is MANDATORY or C++ won't get the expected type.
        if l == k:
            Aij = np.append(Aij, np.float32( 200.)) # Casting value on append is MANDATORY or C++ won't get the expected type.
            Bij = np.append(Bij, np.float32( 33.3)) # Casting value on append is MANDATORY or C++ won't get the expected type.
        if l == k-1 or l == k+1:
            Aij = np.append(Aij, np.float32(-100.)) # Casting value on append is MANDATORY or C++ won't get the expected type.
            Bij = np.append(Bij, np.float32( 16.6)) # Casting value on append is MANDATORY or C++ won't get the expected type.
for k, l, Akl in zip(i, j, Aij):
    print("A[", k, ",", l, "] =", Akl)
for k, l, Bkl in zip(i, j, Bij):
    print("B[", k, ",", l, "] =", Bkl)
A = (n, i, j, Aij) # coo format: dimension, i 1-based indices, j 1-based indices, Aij values.
B = (n, i, j, Bij) # coo format: dimension, i 1-based indices, j 1-based indices, Bij values.

# Get and tune arpack solver.

arpackSlv = pyarpackSlv.float() # Caution: float <=> np.array(..., dtype='float32')
arpackSlv.verbose = 3 # Set to 0 to get a quiet solve.
arpackSlv.debug = 1 # Set to 0 to get a quiet solve.
arpackSlv.nbEV = 2
arpackSlv.nbCV = 2*arpackSlv.nbEV + 1
arpackSlv.mag = 'LM'
arpackSlv.maxIt = 200
arpackSlv.slvTol = 1.e-6
arpackSlv.slvMaxIt = 100
arpackSlv.sigmaReal = 1

# Solve eigen problem.

rc = arpackSlv.solve(A, B)
assert rc == 0, "bad solve"
rc = arpackSlv.checkEigVec(A, B, 1.e-2)
assert rc == 0, "bad checkEigVec"

# Print out results (mode selected, eigen vectors, eigen values, ...).

assert arpackSlv.nbEV == len(arpackSlv.val), "bad result"
print("\nresults:\n")
print("mode selected:", arpackSlv.mode)
print("nb iterations:", arpackSlv.nbIt)
print("Reverse Communication Interface time:", arpackSlv.rciTime, "s")
for val, vec in zip(arpackSlv.val, arpackSlv.vec):
    print("eigen value:", val)
    print("eigen vector:")
    for v in range(n):
        print(vec[v])
