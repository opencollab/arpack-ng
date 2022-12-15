---
title: ARPACK-NG
date: 2022-12-07
update: 2022-12-07
status: stable
docs: none
extpkgs: none
category:
  - LinearAlgebra
tags:
  - matrix
  - ARPACK
---

ARPACK-NG is a collection of Fortran77 subroutines designed to solve large scale
eigenvalue problems.

<a href="https://travis-ci.org/opencollab/arpack-ng"><img src="https://travis-ci.org/opencollab/arpack-ng.svg"/></a><br/>
[![Coverage Status](https://coveralls.io/repos/github/opencollab/arpack-ng/badge.svg?branch=master)](https://coveralls.io/github/opencollab/arpack-ng?branch=master)

## Important Features

- Reverse Communication Interface (RCI).
- Single and Double Precision Real Arithmetic Versions for Symmetric, Non-symmetric, Standard or Generalized Problems.
- Single and Double Precision Complex Arithmetic Versions for Standard or Generalized Problems.
- Routines for Banded Matrices - Standard or Generalized Problems.
- Routines for The Singular Value Decomposition.
- Example driver routines that may be used as templates to implement numerous
- Shift-Invert strategies for all problem types, data types and precision.
- `arpackmm`: utility to test arpack with matrix market files. Note: to run this utility, you need the eigen library (to handle RCI).

## ILP64 support

- Sequential arpack supports [ILP64](https://www.intel.com/content/www/us/en/develop/documentation/onemkl-linux-developer-guide/top/linking-your-application-with-onemkl/linking-in-detail/linking-with-interface-libraries/using-the-ilp64-interface-vs-lp64-interface.html), but, parallel arpack doesn't.
- Reminder: you can NOT mix `ILP64` with `LP64`. If you compile `arpack-ng` with `ILP64` (resp. `LP64`) support, you MUST insure your BLAS/LAPACK is compliant with `ILP64` (resp. `LP64`).
- users: set `INTERFACE64` at configure time.

Note for F77/F90 developers:

- All files which needs `ILP64` support must include `"arpackicb.h"`.
- When coding, use `i_int` (defined in `arpackicb.h`) instead of `c_int`. `i_int` stands for `iso_c_binding int`: it's `#defined` to `c_int` or `c_int64_t` according to the architecture.

Note for C/C++ developers:

- All files which needs `ILP64` support must include `"arpackdef.h"`.
- When coding, use `a_int` (defined in `arpackdef.h`) instead of `int`. Here, `a_int` stands for "architecture int": it's `#defined` to `int` or `int64_t` according to the architecture.

**Example**: to test arpack with sequential `ILP64` MKL assuming you use gnu compilers

```bash
$./bootstrap
export FFLAGS='-DMKL_ILP64 -I/usr/include/mkl'
export FCFLAGS='-DMKL_ILP64 -I/usr/include/mkl'
export LIBS='-Wl,--no-as-needed -L/usr/lib/x86_64-linux-gnu -lmkl_sequential -lmkl_core -lpthread -lm -ldl'
export INTERFACE64=1
./configure --with-blas=mkl_gf_ilp64 --with-lapack=mkl_gf_ilp64
$ make all check
```

To get `ISO_C_BINDING` support:

```bash
./configure --enable-icb
cmake -D ICB=ON
```

- The install will now provide `arpack.h/hpp`, `parpack.h/hpp` and friends.
- Examples of use can be found in `./TESTS` and` ./PARPACK/TESTS/MPI`.

Those who are interested in `ISO_C_BINDING` support can checkout following links 
for more informations.

- <http://fortranwiki.org/fortran/show/ISO_C_BINDING>
- <http://fortranwiki.org/fortran/show/Generating+C+Interfaces>
- <https://www.roguewave.com/sites/rw/files/attachments/StandardizedMixedLanguageProgrammingforCandFortran.pdf>

## Python support

`pyarpack`: python support based on `Boost.Python.Numpy` exposing C++ API.

Check out `./EXAMPLES/PYARPACK/README` for more details.

## About the project

This project started as a joint project between Debian, Octave and Scilab in order to provide a common and maintained version of arpack.
This is now a community project maintained by a few volunteers.
Indeed, no single release has been published by Rice university for the last few years and since many software (Octave, Scilab, R, Matlab...) 
forked it and implemented their own modifications, arpack-ng aims to tackle this by providing a common repository, maintained versions with a testsuite. 
`arpack-ng` is replacing arpack almost everywhere.

## 📁 Directory structure

- You have successfully unbundled ARPACK-NG` and are now in the ARPACK-NG directory that was created for you.

- The directory SRC contains the top level routines including the highest level **reverse communication interface** routines

  - `ssaupd`, `dsaupd`: symmetric single and double precision
  - `snaupd`, `dnaupd`: non-symmetric single and double precision
  - `cnaupd`, `znaupd`: complex non-symmetric single and double precision
  - The headers of these routines contain full documentation of calling sequence and usage.
  - Additional information is given in the `/DOCUMENTS` directory.

- The directory `PARPACK` contains the Parallel ARPACK routines.

- Example driver programs that illustrate all the computational modes, data types and precisions may be found in the EXAMPLES directory. Upon executing the `ls EXAMPLES` command you should see the following directories

  ```bash
  ├── BAND
  ├── COMPLEX
  ├── Makefile.am
  ├── MATRIX_MARKET
  ├── NONSYM
  ├── PYARPACK
  ├── README
  ├── SIMPLE
  ├── SVD
  └── SYM
  ```

  - Example programs for banded, complex, nonsymmetric, symmetric, and singular value decomposition may be found in the directories BAND, COMPLEX, NONSYM, SYM, SVD respectively.
  - Look at the README file for further information.
  - To get started, get into the SIMPLE directory to see example programs that illustrate the use of ARPACK in the simplest modes of operation for the most commonly posed standard eigenvalue problems.

> Example programs for Parallel ARPACK may be found in the directory `PARPACK/EXAMPLES`. Look at the README file for further information.

## Install 🚀

Unlike ARPACK, ARPACK-NG is providing autotools and cmake based build system. In addition, `ARPACK-NG` also provides
`iso_c_binding` support, which enables to call fortran subroutines natively from C or C++.

First, obtain the source code 📥 from github:

```bash
git clone https://github.com/opencollab/arpack-ng.git
cd ./arpack-ng
```

If you prefer the ssh to obtain the source code, then use:

```bash
git clone git@github.com:opencollab/arpack-ng.git
cd ./arpack-ng
```

> Note, It is recommended to install `arpack` at standard location on your system by using your root privilage.


### Using auto-tools

In the source directory, use the following commands to ◉ auto-configure, 📦 build and install the `arpack-ng`.

```bash
sh bootstrap
./configure
make
make check
make install
```

Congratulations 🎉, you have installed `arpack` lib using autotools. 

You can also customize the installation of `arpack` using the autotools. For example,


```bash
LIBSUFFIX="64" ./configure
make all install
```

```bash
INTERFACE64="1" ITF64SUFFIX="ILP64" ./configure
make all install
```

```bash
LIBSUFFIX="64" INTERFACE64="1" ITF64SUFFIX="ILP64" ./configure
make all install
```

To enable `ISO_C_BINDING` support:

```bash
./configure --enable-icb
```

### Using cmake

You can install `ARPACK-NG` by using cmake. If you do not have cmake, then please download the binary from `pip` using:

```bash
python3 -m pip install cmake
which cmake && cmake --version
```

After installing cmake, follow the instruction given below.

Caution: Make sure you are in source directory of ARPACK-NG.

```bash
mkdir build
cd build
cmake -D EXAMPLES=ON -D MPI=ON -D BUILD_SHARED_LIBS=ON ..
make
make install
```

✨ Congratulations, you have installed `arpack` lib using cmake.

> ⚡ Do not miss the `..` in line 3, as it denotes the location of ARPACK-NG source files.

The above-mentioned process will build everything including the examples and 
parallel support using MPI.

You can customize the build by declaring the cmake options during configuration as shown below: 

```bash
cmake -D LIBSUFFIX="64" ..
make all install
```

```bash
cmake -D INTERFACE64=ON -D ITF64SUFFIX="ILP64" ..
make all install
```

```bash
cmake -D LIBSUFFIX="64" -D INTERFACE64=ON -D ITF64SUFFIX="ILP64" ..
make all install
```

To enable `ISO_C_BINDING` support:

```bash
cmake -D ICB=ON
```

## Extras

On mac OS, with GNU compilers, you may need to customize options:

```bash
LIBS="-framework Accelerate" FFLAGS="-ff2c -fno-second-underscore" FCFLAGS="-ff2c -fno-second-underscore" ./configure
```

You can use arpack in your CMake builds by using `ARPACK::ARPACK` target. For example,

```cmake
FIND_PACKAGE(arpack-ng)
ADD_EXECUTABLE(main main.f)
TARGET_INCLUDE_DIRECTORIES(main PUBLIC ARPACK::ARPACK)
TARGET_LINK_LIBRARIES(main ARPACK::ARPACK)
```

To use PARPACK in your Cmake builds, use `PARPACK::PARPACK` target:

```cmake
FIND_PACKAGE(arpack-ng)
ADD_EXECUTABLE(main main.f)
TARGET_INCLUDE_DIRECTORIES(main PUBLIC PARPACK::PARPACK)
TARGET_LINK_LIBRARIES(main PARPACK::PARPACK)
```

Note: Make sure to update `CMAKE_MODULE_PATH` env variable (otheriwse, `find_package` won't find arpack-ng cmake file).

How to use arpack-ng with Intel MKL:

- Let autotools/cmake find MKL for you based on pkg-config files (setting `PKG_CONFIG_PATH`) or cmake options (`BLA_VENDOR=Intel`).
- Refers to the Intel Link Advisor: <https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-link-line-advisor.html>.

## Documentation

Within DOCUMENTS directory there are three files for templates on how to invoke the computational modes of ARPACK.

- ex-sym.doc
- ex-nonsym.doc and
- ex-complex.doc

Also look in the README.MD file for explanations concerning the
other documents.

## Good luck and enjoy 🎊
