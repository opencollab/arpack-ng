#!/usr/bin/env python

import numpy as np
from pyarpack import denseLLT as pyarpackSlv

# Build laplacian.

n = 4
Aij = np.array([], dtype='float64')
for k in range(n):
    for l in range(n):
        if l == k:
            Aij = np.append(Aij, np.float64( 200.)) # Casting value on append is MANDATORY or C++ won't get the expected type.
        elif l == k-1 or l == k+1:
            Aij = np.append(Aij, np.float64(-100.)) # Casting value on append is MANDATORY or C++ won't get the expected type.
        else:
            Aij = np.append(Aij, np.float64(   0.)) # Casting value on append is MANDATORY or C++ won't get the expected type.
for idx, val in enumerate(Aij):
    print("A[", idx, "] =", val)
A = (Aij, False) # raw format: Aij values, row ordered (or not).

# Get and tune arpack solver.

arpackSlv = pyarpackSlv.double() # Caution: double <=> np.array(..., dtype='float64')
arpackSlv.verbose = 3 # Set to 0 to get a quiet solve.
arpackSlv.debug = 1 # Set to 0 to get a quiet solve.
arpackSlv.nbEV = 1
arpackSlv.nbCV = 2*arpackSlv.nbEV + 1
arpackSlv.mag = 'LM'
arpackSlv.maxIt = 200
arpackSlv.slvOffset = 0.
arpackSlv.slvScale = 1.

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
Aij = np.array([], dtype='float32')
Bij = np.array([], dtype='float32')
for k in range(n):
    for l in range(n):
        if l == k:
            Aij = np.append(Aij, np.float32( 200.)) # Casting value on append is MANDATORY or C++ won't get the expected type.
            Bij = np.append(Bij, np.float32( 33.3)) # Casting value on append is MANDATORY or C++ won't get the expected type.
        elif l == k-1 or l == k+1:
            Aij = np.append(Aij, np.float32(-100.)) # Casting value on append is MANDATORY or C++ won't get the expected type.
            Bij = np.append(Bij, np.float32( 16.6)) # Casting value on append is MANDATORY or C++ won't get the expected type.
        else:
            Aij = np.append(Aij, np.float32(   0.)) # Casting value on append is MANDATORY or C++ won't get the expected type.
            Bij = np.append(Bij, np.float32(   0.)) # Casting value on append is MANDATORY or C++ won't get the expected type.
for idx, val in enumerate(Aij):
    print("A[", idx, "] =", val)
for idx, val in enumerate(Bij):
    print("B[", idx, "] =", val)
A = (Aij, False) # raw format: Aij values, row ordered (or not).
B = (Bij,  True) # raw format: Bij values, row ordered (or not).

# Get and tune arpack solver.

arpackSlv = pyarpackSlv.float() # Caution: float <=> np.array(..., dtype='float32')
arpackSlv.verbose = 3 # Set to 0 to get a quiet solve.
arpackSlv.debug = 1 # Set to 0 to get a quiet solve.
arpackSlv.nbEV = 2
arpackSlv.nbCV = 2*arpackSlv.nbEV + 1
arpackSlv.mag = 'LM'
arpackSlv.maxIt = 200
arpackSlv.slvOffset = 0.
arpackSlv.slvScale = 1.

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
