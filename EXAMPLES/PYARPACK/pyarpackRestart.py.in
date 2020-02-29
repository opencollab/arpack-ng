#!/usr/bin/env python

import numpy as np
from pyarpack import sparseBiCGDiag as pyarpackSlv

# Build laplacian.

n = 8
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
arpackSlv.dumpToFile = True # Dump eigen vectors to arpackSolver.*.out files.
arpackSlv.schur = True # Schur vectors and eigenvectors of A are the same if A is a normal matrix.

# Solve eigen problem.

rc = arpackSlv.solve(A)
assert rc == 0, "bad solve"
rc = arpackSlv.checkEigVec(A)
assert rc == 0, "bad checkEigVec"

nbIt1 = arpackSlv.nbIt

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

# Build laplacian (similar-but-different from the previous one).

n = 8
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
            Aij = np.append(Aij, np.float64( 210.)) # Casting value on append is MANDATORY or C++ won't get the expected type.
        if l == k-1 or l == k+1:
            Aij = np.append(Aij, np.float64( -90.)) # Casting value on append is MANDATORY or C++ won't get the expected type.
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
arpackSlv.restartFromFile = True # Restart from eigen vectors found in arpackSolver.*.out files.
arpackSlv.schur = True # Schur vectors and eigenvectors of A are the same if A is a normal matrix.

# Solve eigen problem.

rc = arpackSlv.solve(A)
assert rc == 0, "bad solve"
rc = arpackSlv.checkEigVec(A)
assert rc == 0, "bad checkEigVec"

nbIt2 = arpackSlv.nbIt
assert nbIt2 < nbIt1, "bad restart" # Restart from the first solve to run the second solve for a similar-but-different A.

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
