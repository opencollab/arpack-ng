cmake_minimum_required(VERSION 3.0)

if (NOT DEFINED CMAKE_BUILD_TYPE)
   set (CMAKE_BUILD_TYPE Release CACHE STRING "Build type")
endif ()

project(arpack C Fortran)

set(arpack_ng_MAJOR_VERSION 3)
set(arpack_ng_MINOR_VERSION 9)
set(arpack_ng_PATCH_VERSION 0)
set(arpack_ng_VERSION ${arpack_ng_MAJOR_VERSION}.${arpack_ng_MINOR_VERSION}.${arpack_ng_PATCH_VERSION})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

# set C99 standard
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)

# Adopted from https://github.com/feymark/arpack.git

if (POLICY CMP0042)
  # enable MACOSX_RPATH by default
  cmake_policy (SET CMP0042 NEW)
endif ()

option(BUILD_SHARED_LIBS "Build shared libraries instead of static libraries" ON)
option(MPI "Enable parallel support" OFF)
option(ICB "Enable support for *[ae]upd_c with ISO_C_BINDING" OFF)
option(ICBEXMM "Enable support for matrix market example based on ICB" OFF)
option(PYTHON3 "Enable python3 support" OFF)
set(BOOST_PYTHON_LIBSUFFIX "" CACHE STRING "suffix to add to custom boost python libs")
option(EXAMPLES "Compile ARPACK examples" OFF)
set(LIBSUFFIX ""
    CACHE STRING "suffix to add to ARPACK libraries names")
set(SYMBOLSUFFIX ""
    CACHE STRING "suffix to add to ARPACK, BLAS and LAPACK function names")
option(INTERFACE64 "use the 64-bit integer interface (ILP64) for ARPACK, BLAS and LAPACK")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# We don't want this to run on every build.
option(COVERALLS "Generate coveralls data" OFF)
if (COVERALLS)
    include(Coveralls)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -O0 -fprofile-arcs -ftest-coverage")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O0 -fprofile-arcs -ftest-coverage")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -fprofile-arcs -ftest-coverage")
    # The no space is by design: issue in cmake. See CMP0004.
    set(EXTRA_LDFLAGS "${EXTRA_LDFLAGS}-lgcov")
endif()

function(prefixlist list_name prefix)
    set(${list_name}_TMP)
    foreach(l ${${list_name}})
      list(APPEND ${list_name}_TMP ${prefix}${l} )
    endforeach()
    set(${list_name} "${${list_name}_TMP}" PARENT_SCOPE)
endfunction(prefixlist)

function(examples list_name)
    foreach(l ${${list_name}})
        get_filename_component(lwe ${l} NAME_WE)
        add_executable(${lwe} ${arpackexample_DIR}/${l} ${examples_EXTRA_SRCS})
	target_link_libraries(${lwe} arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
        add_test(NAME "${lwe}_ex" COMMAND ${lwe} WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
    endforeach()
endfunction(examples)

function(pexamples list_name)
    foreach(l ${${list_name}})
        get_filename_component(lwe ${l} NAME_WE)
        add_executable(${lwe} ${parpackexample_DIR}/${l} )
        target_link_libraries(${lwe} parpack arpack MPI::MPI_Fortran)
        add_test(NAME "${lwe}_ex" COMMAND mpirun -n 2 ./${lwe} WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
    endforeach()
endfunction(pexamples)

if (PYTHON3)
    enable_language(C CXX) # Boost requirement.
    set(CMAKE_CXX_STANDARD 14) # Boost requirement.

    find_package(PythonInterp 3 REQUIRED)
    find_package(PythonLibs 3 REQUIRED)
    find_package(Boost COMPONENTS python${BOOST_PYTHON_LIBSUFFIX} numpy${BOOST_PYTHON_LIBSUFFIX} REQUIRED)

    set(ICBEXMM "ON")
endif ()

if (ICBEXMM)
    find_package(Eigen3 3.3 QUIET)
    if (NOT EIGEN3_FOUND) # If not found, piggy-back pkg-config files.
        message(WARNING "CMake didn't find the Eigen3 package. Try to look for pkg-config file...")
        find_package(PkgConfig REQUIRED)
        pkg_check_modules(EIGEN3 REQUIRED eigen3>=3.3)
        set(EIGEN3_INCLUDE_DIR ${EIGEN3_INCLUDE_DIRS})
    endif ()
    find_program (BASH_PROGRAM bash)

    # Look for headers.
    find_path(EIGEN3_ITERATIVE_SOLVER_DIR NAMES IterativeLinearSolvers PATHS ${EIGEN3_INCLUDE_DIR} PATH_SUFFIXES Eigen)
    if (NOT EIGEN3_ITERATIVE_SOLVER_DIR)
        message(FATAL_ERROR "-- Eigen/IterativeLinearSolvers header not found.")
    endif ()
    find_path(EIGEN3_SPARSE_DIR NAMES Sparse PATHS ${EIGEN3_INCLUDE_DIR} PATH_SUFFIXES Eigen)
    if (NOT EIGEN3_SPARSE_DIR)
        message(FATAL_ERROR "-- Eigen/Sparse header not found.")
    endif ()
    find_path(EIGEN3_SPARSELU_SOLVER_DIR NAMES SparseLU PATHS ${EIGEN3_INCLUDE_DIR} PATH_SUFFIXES Eigen)
    if (NOT EIGEN3_SPARSELU_SOLVER_DIR)
        message(FATAL_ERROR "-- Eigen/SparseLU header not found.")
    endif ()
    find_path(EIGEN3_SPARSEQR_SOLVER_DIR NAMES SparseQR PATHS ${EIGEN3_INCLUDE_DIR} PATH_SUFFIXES Eigen)
    if (NOT EIGEN3_SPARSEQR_SOLVER_DIR)
        message(FATAL_ERROR "-- Eigen/SparseQR header not found.")
    endif ()
    find_path(EIGEN3_SPARSECHOLESKY_SOLVER_DIR NAMES SparseCholesky PATHS ${EIGEN3_INCLUDE_DIR} PATH_SUFFIXES Eigen)
    if (NOT EIGEN3_SPARSECHOLESKY_SOLVER_DIR)
        message(FATAL_ERROR "-- Eigen/SparseCholesky header not found.")
    endif ()
    find_path(EIGEN3_DENSE_DIR NAMES Dense PATHS ${EIGEN3_INCLUDE_DIR} PATH_SUFFIXES Eigen)
    if (NOT EIGEN3_DENSE_DIR)
        message(FATAL_ERROR "-- Eigen/Dense header not found.")
    endif ()
    find_path(EIGEN3_DENSELU_DIR NAMES LU PATHS ${EIGEN3_INCLUDE_DIR} PATH_SUFFIXES Eigen)
    if (NOT EIGEN3_DENSELU_DIR)
        message(FATAL_ERROR "-- Eigen/LU header not found.")
    endif ()
    find_path(EIGEN3_DENSEQR_DIR NAMES QR PATHS ${EIGEN3_INCLUDE_DIR} PATH_SUFFIXES Eigen)
    if (NOT EIGEN3_DENSEQR_DIR)
        message(FATAL_ERROR "-- Eigen/QR header not found.")
    endif ()
    find_path(EIGEN3_DENSECHOLESKY_DIR NAMES Cholesky PATHS ${EIGEN3_INCLUDE_DIR} PATH_SUFFIXES Eigen)
    if (NOT EIGEN3_DENSECHOLESKY_DIR)
        message(FATAL_ERROR "-- Eigen/Cholesky header not found.")
    endif ()

    set(ICB "ON")
endif ()

# Enable language(s) before any find_package (in particular before MPI find_package).
if (ICB)
    enable_language(C CXX) # For testing binding with c/c++.

    set(CMAKE_CXX_STANDARD 11) # OK, since cmake-3.1 only.

    file(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/PROG_ICB.f90
         "
          PROGRAM PROG_ICB
          USE iso_c_binding
          IMPLICIT NONE
          INTEGER(C_INT) :: a
          a = 1
          END PROGRAM PROG_ICB
         ")
    try_compile(COMPILE_ICB ${CMAKE_BINARY_DIR} ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/PROG_ICB.f90)
    if(NOT ${COMPILE_ICB})
      message(FATAL_ERROR "-- Fortran compiler does not support iso_c_binding.")
    else()
      message("-- Fortran compiler does support iso_c_binding.")
    endif()
else()
    # ICB saves you from old-fashion-boring-cumbersome-fortran/C crap... For arpack symbols (only).

    if (SYMBOLSUFFIX)
        if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
          set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -ffixed-line-length-none")
        elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
          set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpp -extend-source")
        else ()
          message(WARNING "build script does not know how to enable your Fortran compiler's preprocessor and support for lines longer than 72 characters: set them manually via FFLAGS.")
        endif ()

        set(SCALARFUNS naitr napps naup2 naupd nconv neigh neupd ngets statn saitr sapps saup2 saupd sconv seigt seupd sgets stats getv0 sortc sortr sesrt stqrb)
        set(COMPLEXFUNS, naitr napps naup2 naupd neigh neupd ngets statn getv0 sortc)

        foreach (f IN LISTS SCALARFUNS)
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Ds${f}=s${f}${SYMBOLSUFFIX} -Dd${f}=d${f}${SYMBOLSUFFIX}")
        endforeach ()

        foreach (f IN LISTS COMPLEXFUNS)
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Dc${f}=c${f}${SYMBOLSUFFIX} -Dz${f}=z${f}${SYMBOLSUFFIX}")
        endforeach ()

        set(CFUNS snaupd sneupd dnaupd dneupd)
        foreach (f IN LISTS CFUNS)
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D${f}=${f}${SYMBOLSUFFIX}")
            list(APPEND CFUNS_SUFFIXED ${f}${SYMBOLSUFFIX})
        endforeach ()

        include(FortranCInterface)
        FortranCInterface_HEADER(FCMangle.h SYMBOLS ${CFUNS_SUFFIXED})

        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DINCLUDE_FCMANGLE")

        FortranCInterface_VERIFY()
    endif ()
endif ()

# Find MPI

if (MPI)
    if (NOT TARGET MPI::MPI_Fortran) # Search only if not already found by upper CMakeLists.txt
        include(FindMPI)
        find_package(MPI REQUIRED COMPONENTS Fortran)

        # MPI::MPI_* target was already created at this point by FindMPI.cmake if cmake version >= 3.9
        if (NOT TARGET MPI::MPI_Fortran) # Create target "at hand" to ensure compatibility if cmake version < 3.9
            add_library(MPI::MPI_Fortran INTERFACE IMPORTED)
            set_target_properties(MPI::MPI_Fortran PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${MPI_Fortran_INCLUDE_DIRS}")
            set_target_properties(MPI::MPI_Fortran PROPERTIES INTERFACE_LINK_LIBRARIES      "${MPI_Fortran_LIBRARIES}")
        endif()
    endif()
    get_target_property(MPI_Fortran_INCLUDE_DIRS MPI::MPI_Fortran INTERFACE_INCLUDE_DIRECTORIES) # Get variables from target (*.pc/cmake, msg).
    get_target_property(MPI_Fortran_LIBRARIES    MPI::MPI_Fortran INTERFACE_LINK_LIBRARIES)      # Get variables from target (*.pc/cmake, msg).

    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAG}")

    if(ICB)
        if (NOT TARGET MPI::MPI_C) # Search only if not already found by upper CMakeLists.txt
            include(FindMPI)
            find_package(MPI REQUIRED COMPONENTS C)

            if (NOT TARGET MPI::MPI_C) # Create target "at hand" to ensure compatibility if cmake version < 3.9
                add_library(MPI::MPI_C INTERFACE IMPORTED)
                set_target_properties(MPI::MPI_C PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${MPI_C_INCLUDE_DIRS}")
                set_target_properties(MPI::MPI_C PROPERTIES INTERFACE_LINK_LIBRARIES      "${MPI_C_LIBRARIES}")
            endif()
        endif()
        get_target_property(MPI_C_INCLUDE_DIRS   MPI::MPI_C       INTERFACE_INCLUDE_DIRECTORIES) # Get variables from target (*.pc/cmake, msg).
        get_target_property(MPI_C_LIBRARIES      MPI::MPI_C       INTERFACE_LINK_LIBRARIES)      # Get variables from target (*.pc/cmake, msg).

        if (NOT TARGET MPI::MPI_CXX) # Search only if not already found by upper CMakeLists.txt
            include(FindMPI)
            find_package(MPI REQUIRED COMPONENTS CXX)

            if (NOT TARGET MPI::MPI_CXX) # Create target "at hand" to ensure compatibility if cmake version < 3.9
                add_library(MPI::MPI_CXX INTERFACE IMPORTED)
                set_target_properties(MPI::MPI_CXX PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${MPI_CXX_INCLUDE_DIRS}")
                set_target_properties(MPI::MPI_CXX PROPERTIES INTERFACE_LINK_LIBRARIES      "${MPI_CXX_LIBRARIES}")
            endif()
        endif()
        get_target_property(MPI_CXX_INCLUDE_DIRS MPI::MPI_CXX     INTERFACE_INCLUDE_DIRECTORIES) # Get variables from target (*.pc/cmake, msg).
        get_target_property(MPI_CXX_LIBRARIES    MPI::MPI_CXX     INTERFACE_LINK_LIBRARIES)      # Get variables from target (*.pc/cmake, msg).

        include(CheckSymbolExists)
        check_symbol_exists(MPI_Comm_c2f "${MPI_C_INCLUDE_DIRS}/mpi.h" MPI_Comm_c2f_FOUND)
        if(NOT ${MPI_Comm_c2f_FOUND})
          message(FATAL_ERROR "symbol MPI_Comm_c2f does not exist")
        endif()
    endif()
endif()

# Find BLAS

if (NOT TARGET BLAS::BLAS) # Search only if not already found by upper CMakeLists.txt
    find_package(BLAS REQUIRED)

    # BLAS::BLAS target was already created at this point by FindBLAS.cmake if cmake version >= 3.18
    if (NOT TARGET BLAS::BLAS) # Create target "at hand" to ensure compatibility if cmake version < 3.18
        add_library(BLAS::BLAS INTERFACE IMPORTED)
        set_target_properties(BLAS::BLAS PROPERTIES INTERFACE_LINK_LIBRARIES "${BLAS_LIBRARIES}")
    endif()
endif()
get_target_property(BLAS_LIBRARIES BLAS::BLAS INTERFACE_LINK_LIBRARIES) # Get variables from target (*.pc/cmake, msg).

# Find LAPACK

if (NOT TARGET LAPACK::LAPACK) # Search only if not already found by upper CMakeLists.txt
    find_package(LAPACK REQUIRED)

    # LAPACK::LAPACK target was already created at this point by FindLAPACK.cmake if cmake version >= 3.18
    if (NOT TARGET LAPACK::LAPACK) # Create target "at hand" to ensure compatibility if cmake version < 3.18
        add_library(LAPACK::LAPACK INTERFACE IMPORTED)
        set_target_properties(LAPACK::LAPACK PROPERTIES INTERFACE_LINK_LIBRARIES "${LAPACK_LIBRARIES}")
    endif()
endif()
get_target_property(LAPACK_LIBRARIES LAPACK::LAPACK INTERFACE_LINK_LIBRARIES) # Get variables from target (*.pc/cmake, msg).

# As BLAS/LAPACK does not provide ICB, we may have to deal with symbols the old-fashion-boring-cumbersome-fortran/C way...

if (SYMBOLSUFFIX)
    if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -ffixed-line-length-none")
    elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpp -extend-source")
    else ()
      message(WARNING "build script does not know how to enable your Fortran compiler's preprocessor and support for lines longer than 72 characters: set them manually via FFLAGS.")
    endif ()

    set(BLASFUNS1 axpy copy gemv geqr2 lacpy lae2 lahqr lanhs larnv lartg lascl laset lasrt scal trevc trmm trsen gbmv gbtrf gbtrs gttrf gttrs pttrf pttrs)
    set(BLASFUNS2 dot ger labad laev2 lamch lanst lanv2 lapy2 larf larfg lasr nrm2 orm2r rot steqr swap)
    set(BLASFUNS3 dotc geru unm2r)
    set(BLASFUNS4 COPY LABAD LAMCH LANHS LANV2 LARFG ROT GEMV)
    set(BLASFUNS5 scnrm2 dznrm2 csscal zdscal)

    foreach (f IN LISTS BLASFUNS1 BLASFUNS2)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Ds${f}=s${f}${SYMBOLSUFFIX} -Dd${f}=d${f}${SYMBOLSUFFIX}")
    endforeach ()

    foreach (f IN LISTS BLASFUNS1 BLASFUNS3)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Dc${f}=c${f}${SYMBOLSUFFIX} -Dz${f}=z${f}${SYMBOLSUFFIX}")
    endforeach ()

    foreach (f IN LISTS BLASFUNS4)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DS${f}=S${f}${SYMBOLSUFFIX} -DD${f}=D${f}${SYMBOLSUFFIX}")
    endforeach ()

    foreach (f IN LISTS BLASFUNS5)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -D${f}=${f}${SYMBOLSUFFIX}")
    endforeach ()

    set(CFUNS sgemm cheev)
    foreach (f IN LISTS CFUNS)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D${f}=${f}${SYMBOLSUFFIX}")
        list(APPEND CFUNS_SUFFIXED ${f}${SYMBOLSUFFIX})
    endforeach ()

    include(FortranCInterface)
    FortranCInterface_HEADER(FCMangle.h SYMBOLS ${CFUNS_SUFFIXED})

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DINCLUDE_FCMANGLE")

    FortranCInterface_VERIFY()
endif ()

if (MPI)
    set(parpackutil_DIR ${arpack_SOURCE_DIR}/PARPACK/UTIL/)
    set(parpacksrc_DIR ${arpack_SOURCE_DIR}/PARPACK/SRC/)
endif()

file(GLOB arpacksrc_STAT_SRCS ${arpack_SOURCE_DIR}/dbgini.f ${arpack_SOURCE_DIR}/staini.f ${arpack_SOURCE_DIR}/SRC/*.f)

set(arpacksrc_ICB "")
set(parpacksrc_ICB "")
if(ICB)
  file(GLOB arpacksrc_ICB ${arpack_SOURCE_DIR}/SRC/icba*.F90 ${arpack_SOURCE_DIR}/ICB/debug_icb.F90 ${arpack_SOURCE_DIR}/ICB/stat_icb.F90)
  file(GLOB parpacksrc_ICB ${arpack_SOURCE_DIR}/PARPACK/SRC/MPI/icbp*.F90 ${arpack_SOURCE_DIR}/ICB/debug_icb.F90 ${arpack_SOURCE_DIR}/ICB/stat_icb.F90)
endif()

set(arpackutil_STAT_SRCS
    ${arpack_SOURCE_DIR}/UTIL/icnteq.f
    ${arpack_SOURCE_DIR}/UTIL/icopy.f
    ${arpack_SOURCE_DIR}/UTIL/iset.f
    ${arpack_SOURCE_DIR}/UTIL/iswap.f
    ${arpack_SOURCE_DIR}/UTIL/ivout.f
    ${arpack_SOURCE_DIR}/UTIL/second_NONE.f
    ${arpack_SOURCE_DIR}/UTIL/svout.f
    ${arpack_SOURCE_DIR}/UTIL/smout.f
    ${arpack_SOURCE_DIR}/UTIL/dvout.f
    ${arpack_SOURCE_DIR}/UTIL/dmout.f
    ${arpack_SOURCE_DIR}/UTIL/cvout.f
    ${arpack_SOURCE_DIR}/UTIL/cmout.f
    ${arpack_SOURCE_DIR}/UTIL/zvout.f
    ${arpack_SOURCE_DIR}/UTIL/zmout.f )


if (MPI)
    file(GLOB parpackutil_STAT_SRCS ${arpack_SOURCE_DIR}/PARPACK/UTIL/MPI/*.f)
    file(GLOB parpacksrc_STAT_SRCS ${arpack_SOURCE_DIR}/dbgini.f ${arpack_SOURCE_DIR}/staini.f ${arpack_SOURCE_DIR}/PARPACK/SRC/MPI/*.f)
endif()

# use -DBUILD_SHARED_LIBS=ON|OFF to control static/shared
add_library(arpack ${arpackutil_STAT_SRCS} ${arpacksrc_STAT_SRCS} ${arpacksrc_ICB})

target_link_libraries(arpack BLAS::BLAS)
target_link_libraries(arpack LAPACK::LAPACK ${EXTRA_LDFLAGS})
set_target_properties(arpack PROPERTIES OUTPUT_NAME arpack${LIBSUFFIX})
set_target_properties(arpack PROPERTIES VERSION 2.1.0)
set_target_properties(arpack PROPERTIES SOVERSION 2)

if (MPI)
    # use -DBUILD_SHARED_LIBS=ON|OFF to control static/shared
    add_library(parpack ${parpacksrc_STAT_SRCS} ${parpackutil_STAT_SRCS} ${parpacksrc_ICB})

    target_link_libraries(parpack MPI::MPI_Fortran)
    target_link_libraries(parpack arpack)
    set_target_properties(parpack PROPERTIES OUTPUT_NAME parpack${LIBSUFFIX})
    set_target_properties(parpack PROPERTIES VERSION 2.1.0)
    set_target_properties(parpack PROPERTIES SOVERSION 2)
endif ()

set(PYINT "int32" CACHE STRING "int type to be used in python scripts") # PYINT : int used with python.
if (INTERFACE64)
    set(PYINT "int64" CACHE STRING "int type to be used in python scripts") # PYINT : int used with python.
    set(INTERFACE64 1)
    if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-integer-8")
    elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -i8")
    else ()
      message(WARNING "build script does not know how to make your Fortran compiler use 64-bit integers: set it manually via FFLAGS.")
    endif ()
else ()
    set(INTERFACE64 0)
endif ()

include_directories("${PROJECT_BINARY_DIR}") # Find arpackdef.h and arpackicb.h
set(CMAKE_INCLUDE_CURRENT_DIR ON)

############################
# EXAMPLES
############################
if (EXAMPLES)
    ############################
    # EXAMPLES/BAND
    ############################
    set(arpackexample_DIR ${arpack_SOURCE_DIR}/EXAMPLES/BAND/)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/EXAMPLES/BAND/)

    set(examples_EXTRA_SRCS ${arpackexample_DIR}/cnband.f)
    set(examples_STAT_SRCS
        cnbdr1.f
        cnbdr2.f
        cnbdr3.f
        cnbdr4.f)
    examples(examples_STAT_SRCS)

    set(examples_EXTRA_SRCS ${arpackexample_DIR}/dnband.f)
    set(examples_STAT_SRCS
        dnbdr1.f
        dnbdr2.f
        dnbdr3.f
        dnbdr4.f
        dnbdr5.f
        dnbdr6.f)
    examples(examples_STAT_SRCS)

    set(examples_EXTRA_SRCS ${arpackexample_DIR}/dsband.f)
    set(examples_STAT_SRCS
        dsbdr1.f
        dsbdr2.f
        dsbdr3.f
        dsbdr4.f
        dsbdr5.f
        dsbdr6.f)
    examples(examples_STAT_SRCS)

    set(examples_EXTRA_SRCS ${arpackexample_DIR}/snband.f)
    set(examples_STAT_SRCS
        snbdr1.f
        snbdr2.f
        snbdr3.f
        snbdr4.f
        snbdr5.f
        snbdr6.f)
    examples(examples_STAT_SRCS)

    set(examples_EXTRA_SRCS ${arpackexample_DIR}/ssband.f)
    set(examples_STAT_SRCS
        ssbdr1.f
        ssbdr2.f
        ssbdr3.f
        ssbdr4.f
        ssbdr5.f
        ssbdr6.f)
    examples(examples_STAT_SRCS)

    set(examples_EXTRA_SRCS ${arpackexample_DIR}/znband.f)
    set(examples_STAT_SRCS
        znbdr1.f
        znbdr2.f
        znbdr3.f
        znbdr4.f)
    examples(examples_STAT_SRCS)
    set(examples_EXTRA_SRCS)

    ############################
    # EXAMPLES/COMPLEX
    ############################
    set(arpackexample_DIR ${arpack_SOURCE_DIR}/EXAMPLES/COMPLEX/)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/EXAMPLES/COMPLEX/)

    set(examples_STAT_SRCS
        cndrv1.f
        cndrv2.f
        cndrv3.f
        cndrv4.f
        zndrv1.f
        zndrv2.f
        zndrv3.f
        zndrv4.f)
    examples(examples_STAT_SRCS)

    ############################
    # examples/nonsym
    ############################
    set(arpackexample_DIR ${arpack_SOURCE_DIR}/EXAMPLES/NONSYM/)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/EXAMPLES/NONSYM/)

    set(examples_STAT_SRCS
        dndrv1.f
        dndrv2.f
        dndrv3.f
        dndrv4.f
        dndrv5.f
        dndrv6.f
        sndrv1.f
        sndrv2.f
        sndrv3.f
        sndrv4.f
        sndrv5.f
        sndrv6.f)
    examples(examples_STAT_SRCS)

    ############################
    # examples/SIMPLE
    ############################
    set(arpackexample_DIR ${arpack_SOURCE_DIR}/EXAMPLES/SIMPLE/)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/EXAMPLES/SIMPLE/)

    set(examples_STAT_SRCS
        cnsimp.f
        dnsimp.f
        dssimp.f
        snsimp.f
        sssimp.f
        znsimp.f)

    examples(examples_STAT_SRCS)

    ############################
    # examples/svd
    ############################
    set(arpackexample_DIR ${arpack_SOURCE_DIR}/EXAMPLES/SVD/)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/EXAMPLES/SVD/)

    set(examples_STAT_SRCS
        ssvd.f)

    examples(examples_STAT_SRCS)

    ############################
    # examples/sym
    ############################
    set(arpackexample_DIR ${arpack_SOURCE_DIR}/EXAMPLES/SYM/)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/EXAMPLES/SYM/)

    set(examples_STAT_SRCS
        dsdrv1.f
        dsdrv2.f
        dsdrv3.f
        dsdrv4.f
        dsdrv5.f
        dsdrv6.f
        ssdrv1.f
        ssdrv2.f
        ssdrv3.f
        ssdrv4.f
        ssdrv5.f
        ssdrv6.f)

    examples(examples_STAT_SRCS)

    ############################
    # parpack/examples
    ############################
    if (MPI)
        set(parpackexample_DIR ${arpack_SOURCE_DIR}/PARPACK/EXAMPLES/MPI/)
        set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/PARPACK/EXAMPLES/MPI/)

        set(pexamples_STAT_SRCS
            pcndrv1.f
            pdndrv1.f
            pdndrv3.f
            pdsdrv1.f
            psndrv1.f
            psndrv3.f
            pssdrv1.f
            pzndrv1.f)

        pexamples(pexamples_STAT_SRCS)
    endif()
endif()


############################
# TEST
############################

enable_testing()

set(CMAKE_CTEST_COMMAND ctest -V)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/TESTS)

add_executable(dnsimp_test TESTS/dnsimp.f TESTS/mmio.f TESTS/debug.h)
set_target_properties( dnsimp_test PROPERTIES OUTPUT_NAME  dnsimp )
target_link_libraries(dnsimp_test arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
add_custom_command(TARGET dnsimp_test POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/TESTS/testA.mtx testA.mtx
)
add_test(dnsimp_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/dnsimp)

if (ICB)
    add_executable(bug_1315_single TESTS/bug_1315_single.c)
    target_include_directories(bug_1315_single PUBLIC ${PROJECT_SOURCE_DIR}/ICB) # Get arpack.h
    target_link_libraries(bug_1315_single arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
    add_test(bug_1315_single_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/bug_1315_single)

    add_executable(bug_1315_double TESTS/bug_1315_double.c)
    target_include_directories(bug_1315_double PUBLIC ${PROJECT_SOURCE_DIR}/ICB) # Get arpack.h
    target_link_libraries(bug_1315_double arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
    add_test(bug_1315_double_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/bug_1315_double)
endif()

add_executable(bug_1323 TESTS/bug_1323.f)
target_link_libraries(bug_1323 arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
add_test(bug_1323_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/bug_1323)

add_executable(bug_58_double TESTS/bug_58_double.f)
target_link_libraries(bug_58_double arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
add_test(bug_58_double_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/bug_58_double)

add_executable(bug_79_double_complex TESTS/bug_79_double_complex.f)
target_link_libraries(bug_79_double_complex arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
add_test(bug_79_double_complex_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/bug_79_double_complex)

add_executable(bug_142 TESTS/bug_142.f)
target_link_libraries(bug_142 arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
add_test(bug_142_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/bug_142)

add_executable(bug_142_gen TESTS/bug_142_gen.f)
target_link_libraries(bug_142_gen arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
add_test(bug_142_gen_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/bug_142_gen)

if(MPI)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/PARPACK/TESTS/MPI)

  add_executable(issue46 PARPACK/TESTS/MPI/issue46.f)
  target_link_libraries(issue46 parpack arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
  add_test(issue46_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/issue46)
endif()

if(ICB)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/TESTS)

  add_executable(icb_arpack_c TESTS/icb_arpack_c.c)
  target_include_directories(icb_arpack_c PUBLIC ${PROJECT_SOURCE_DIR}/ICB) # Get arpack.h
  target_link_libraries(icb_arpack_c arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
  add_test(icb_arpack_c_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/icb_arpack_c)

  add_executable(icb_arpack_cpp TESTS/icb_arpack_cpp.cpp)
  target_include_directories(icb_arpack_cpp PUBLIC ${PROJECT_SOURCE_DIR}/ICB) # Get arpack.hpp
  target_link_libraries(icb_arpack_cpp arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
  add_test(icb_arpack_cpp_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/icb_arpack_cpp)

  if (ICBEXMM)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/EXAMPLES/MATRIX_MARKET)

    add_executable(arpackmm EXAMPLES/MATRIX_MARKET/arpackmm.cpp)
    target_include_directories(arpackmm PUBLIC ${PROJECT_SOURCE_DIR}/ICB ${EIGEN3_INCLUDE_DIR}) # Get arpack.h + eigen
    target_link_libraries(arpackmm arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
    configure_file(EXAMPLES/MATRIX_MARKET/As.mtx ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/As.mtx)
    configure_file(EXAMPLES/MATRIX_MARKET/An.mtx ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/An.mtx)
    configure_file(EXAMPLES/MATRIX_MARKET/Az.mtx ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Az.mtx)
    configure_file(EXAMPLES/MATRIX_MARKET/B.mtx  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/B.mtx)
    configure_file(EXAMPLES/MATRIX_MARKET/Bz.mtx ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Bz.mtx)
    configure_file(EXAMPLES/MATRIX_MARKET/arpackmm.sh ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/arpackmm.sh)
    add_test(NAME arpackmm_tst WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} COMMAND ${BASH_PROGRAM} arpackmm.sh)
  endif()

  if (PYTHON3)
    python_add_module(pyarpack ${arpackutil_STAT_SRCS} ${arpacksrc_STAT_SRCS} ${arpacksrc_ICB} ${PROJECT_SOURCE_DIR}/EXAMPLES/PYARPACK/pyarpack.cpp)
    target_compile_definitions(pyarpack PRIVATE PY_MAJOR_VERSION="3")
    set(pyarpack_HDR ${PROJECT_SOURCE_DIR}/ICB ${PROJECT_SOURCE_DIR}/EXAMPLES/MATRIX_MARKET ${PROJECT_SOURCE_DIR}/EXAMPLES/PYARPACK)
    target_include_directories(pyarpack PUBLIC ${pyarpack_HDR} ${EIGEN3_INCLUDE_DIR} ${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS})
    target_link_libraries(pyarpack BLAS::BLAS LAPACK::LAPACK ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
    install(TARGETS pyarpack
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/pyarpack
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/pyarpack)
    configure_file("${PROJECT_SOURCE_DIR}/EXAMPLES/PYARPACK/pyarpackSparseBiCGDiag.py.in" "${CMAKE_BINARY_DIR}/pyarpackSparseBiCGDiag.py" @ONLY)
    add_test(NAME pyarpackSparseBiCGDiag_tst COMMAND ${PYTHON_EXECUTABLE} pyarpackSparseBiCGDiag.py)
    set_tests_properties(pyarpackSparseBiCGDiag_tst PROPERTIES ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/lib:$ENV{PYTHONPATH})
    configure_file("${PROJECT_SOURCE_DIR}/EXAMPLES/PYARPACK/pyarpackSparseBiCGILU.py.in"  "${CMAKE_BINARY_DIR}/pyarpackSparseBiCGILU.py" @ONLY)
    add_test(NAME pyarpackSparseBiCGILU_tst COMMAND ${PYTHON_EXECUTABLE} pyarpackSparseBiCGILU.py)
    set_tests_properties(pyarpackSparseBiCGILU_tst PROPERTIES ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/lib:$ENV{PYTHONPATH})
    configure_file("${PROJECT_SOURCE_DIR}/EXAMPLES/PYARPACK/pyarpackSparseCGDiag.py.in" "${CMAKE_BINARY_DIR}/pyarpackSparseCGDiag.py" @ONLY)
    add_test(NAME pyarpackSparseCGDiag_tst COMMAND ${PYTHON_EXECUTABLE} pyarpackSparseCGDiag.py)
    set_tests_properties(pyarpackSparseCGDiag_tst PROPERTIES ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/lib:$ENV{PYTHONPATH})
    configure_file("${PROJECT_SOURCE_DIR}/EXAMPLES/PYARPACK/pyarpackSparseCGILU.py.in" "${CMAKE_BINARY_DIR}/pyarpackSparseCGILU.py" @ONLY)
    add_test(NAME pyarpackSparseCGILU_tst COMMAND ${PYTHON_EXECUTABLE} pyarpackSparseCGILU.py)
    set_tests_properties(pyarpackSparseCGILU_tst PROPERTIES ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/lib:$ENV{PYTHONPATH})
    configure_file("${PROJECT_SOURCE_DIR}/EXAMPLES/PYARPACK/pyarpackSparseLLT.py.in" "${CMAKE_BINARY_DIR}/pyarpackSparseLLT.py" @ONLY)
    add_test(NAME pyarpackSparseLLT_tst COMMAND ${PYTHON_EXECUTABLE} pyarpackSparseLLT.py)
    set_tests_properties(pyarpackSparseLLT_tst PROPERTIES ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/lib:$ENV{PYTHONPATH})
    configure_file("${PROJECT_SOURCE_DIR}/EXAMPLES/PYARPACK/pyarpackSparseLDLT.py.in" "${CMAKE_BINARY_DIR}/pyarpackSparseLDLT.py" @ONLY)
    add_test(NAME pyarpackSparseLDLT_tst COMMAND ${PYTHON_EXECUTABLE} pyarpackSparseLDLT.py)
    set_tests_properties(pyarpackSparseLDLT_tst PROPERTIES ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/lib:$ENV{PYTHONPATH})
    configure_file("${PROJECT_SOURCE_DIR}/EXAMPLES/PYARPACK/pyarpackSparseLU.py.in" "${CMAKE_BINARY_DIR}/pyarpackSparseLU.py" @ONLY)
    add_test(NAME pyarpackSparseLU_tst COMMAND ${PYTHON_EXECUTABLE} pyarpackSparseLU.py)
    set_tests_properties(pyarpackSparseLU_tst PROPERTIES ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/lib:$ENV{PYTHONPATH})
    configure_file("${PROJECT_SOURCE_DIR}/EXAMPLES/PYARPACK/pyarpackSparseQR.py.in" "${CMAKE_BINARY_DIR}/pyarpackSparseQR.py" @ONLY)
    add_test(NAME pyarpackSparseQR_tst COMMAND ${PYTHON_EXECUTABLE} pyarpackSparseQR.py)
    set_tests_properties(pyarpackSparseQR_tst PROPERTIES ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/lib:$ENV{PYTHONPATH})
    configure_file("${PROJECT_SOURCE_DIR}/EXAMPLES/PYARPACK/pyarpackDenseLLT.py.in" "${CMAKE_BINARY_DIR}/pyarpackDenseLLT.py" @ONLY)
    add_test(NAME pyarpackDenseLLT_tst COMMAND ${PYTHON_EXECUTABLE} pyarpackDenseLLT.py)
    set_tests_properties(pyarpackDenseLLT_tst PROPERTIES ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/lib:$ENV{PYTHONPATH})
    configure_file("${PROJECT_SOURCE_DIR}/EXAMPLES/PYARPACK/pyarpackDenseLDLT.py.in" "${CMAKE_BINARY_DIR}/pyarpackDenseLDLT.py" @ONLY)
    add_test(NAME pyarpackDenseLDLT_tst COMMAND ${PYTHON_EXECUTABLE} pyarpackDenseLDLT.py)
    set_tests_properties(pyarpackDenseLDLT_tst PROPERTIES ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/lib:$ENV{PYTHONPATH})
    configure_file("${PROJECT_SOURCE_DIR}/EXAMPLES/PYARPACK/pyarpackDenseLURR.py.in" "${CMAKE_BINARY_DIR}/pyarpackDenseLURR.py" @ONLY)
    add_test(NAME pyarpackDenseLURR_tst COMMAND ${PYTHON_EXECUTABLE} pyarpackDenseLURR.py)
    set_tests_properties(pyarpackDenseLURR_tst PROPERTIES ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/lib:$ENV{PYTHONPATH})
    configure_file("${PROJECT_SOURCE_DIR}/EXAMPLES/PYARPACK/pyarpackDenseQRRR.py.in" "${CMAKE_BINARY_DIR}/pyarpackDenseQRRR.py" @ONLY)
    add_test(NAME pyarpackDenseQRRR_tst COMMAND ${PYTHON_EXECUTABLE} pyarpackDenseQRRR.py)
    set_tests_properties(pyarpackDenseQRRR_tst PROPERTIES ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/lib:$ENV{PYTHONPATH})
    configure_file("${PROJECT_SOURCE_DIR}/EXAMPLES/PYARPACK/pyarpackDenseLUPP.py.in" "${CMAKE_BINARY_DIR}/pyarpackDenseLUPP.py" @ONLY)
    add_test(NAME pyarpackDenseLUPP_tst COMMAND ${PYTHON_EXECUTABLE} pyarpackDenseLUPP.py)
    set_tests_properties(pyarpackDenseLUPP_tst PROPERTIES ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/lib:$ENV{PYTHONPATH})
    configure_file("${PROJECT_SOURCE_DIR}/EXAMPLES/PYARPACK/pyarpackDenseQRPP.py.in" "${CMAKE_BINARY_DIR}/pyarpackDenseQRPP.py" @ONLY)
    add_test(NAME pyarpackDenseQRPP_tst COMMAND ${PYTHON_EXECUTABLE} pyarpackDenseQRPP.py)
    set_tests_properties(pyarpackDenseQRPP_tst PROPERTIES ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/lib:$ENV{PYTHONPATH})
    configure_file("${PROJECT_SOURCE_DIR}/EXAMPLES/PYARPACK/pyarpackRestart.py.in" "${CMAKE_BINARY_DIR}/pyarpackRestart.py" @ONLY)
    add_test(NAME pyarpackRestart_tst COMMAND ${PYTHON_EXECUTABLE} pyarpackRestart.py)
    set_tests_properties(pyarpackRestart_tst PROPERTIES ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/lib:$ENV{PYTHONPATH})
  endif()

  if (MPI)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/PARPACK/TESTS/MPI)

    add_executable(icb_parpack_c PARPACK/TESTS/MPI/icb_parpack_c.c)
    target_include_directories(icb_parpack_c PUBLIC ${PROJECT_SOURCE_DIR}/ICB MPI::MPI_C) # Get parpack.h mpi.h
    target_link_libraries(icb_parpack_c parpack arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS} MPI::MPI_C)
    add_test(icb_parpack_c_tst mpirun -n 2 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/icb_parpack_c)

    add_executable(icb_parpack_cpp PARPACK/TESTS/MPI/icb_parpack_cpp.cpp)
    target_include_directories(icb_parpack_cpp PUBLIC ${PROJECT_SOURCE_DIR}/ICB MPI::MPI_CXX) # Get parpack.hpp mpi.h
    target_link_libraries(icb_parpack_cpp parpack arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS} MPI::MPI_CXX)
    add_test(icb_parpack_cpp_tst mpirun -n 2 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/icb_parpack_cpp)
  endif()
endif()

############################
# install
############################
# 'make install' to the correct location
include(GNUInstallDirs)

# Convert variable names to those expected by the .pc file.
set(prefix ${CMAKE_INSTALL_PREFIX})
set(exec_prefix \${prefix})
set(libdir ${CMAKE_INSTALL_FULL_LIBDIR})
set(includedir ${CMAKE_INSTALL_FULL_INCLUDEDIR})
set(PACKAGE_NAME ${PROJECT_NAME})
set(PACKAGE_VERSION ${arpack_ng_VERSION})
set(PACKAGE_URL "https://github.com/opencollab/arpack-ng/")

# Convert (LAPACK|BLAS)_LIBRARIES to -l flags, store in ARPACK_PC_LIBS_PRIVATE.
# LAPACK_LIBRARIES contains transitive deps, no need to parse BLAS_LIBRARIES.
set(ARPACK_PC_LIBS_PRIVATE)
foreach(lib ${LAPACK_LIBRARIES})
    get_filename_component(libname ${lib} NAME)
    string(REGEX REPLACE "^lib([^.]+).*$" "-l\\1" libname ${libname})
    list(APPEND ARPACK_PC_LIBS_PRIVATE "${libname}")
endforeach()
string(REPLACE ";" " " ARPACK_PC_LIBS_PRIVATE "${ARPACK_PC_LIBS_PRIVATE}")

set(PARPACK_PC_LIBS_PRIVATE)
foreach(lib ${LAPACK_LIBRARIES} ${MPI_Fortran_LIBRARIES})
    get_filename_component(libname ${lib} NAME)
    string(REGEX REPLACE "^lib([^.]+).*$" "-l\\1" libname ${libname})
    list(APPEND PARPACK_PC_LIBS_PRIVATE "${libname}")
endforeach()
string(REPLACE ";" " " PARPACK_PC_LIBS_PRIVATE "${PARPACK_PC_LIBS_PRIVATE}")

configure_file(SRC/arpack.pc.in "${PROJECT_BINARY_DIR}/SRC/arpack${LIBSUFFIX}.pc" @ONLY)
configure_file(PARPACK/SRC/MPI/parpack.pc.in "${PROJECT_BINARY_DIR}/PARPACK/SRC/MPI/parpack${LIBSUFFIX}.pc" @ONLY)
configure_file(EXAMPLES/MATRIX_MARKET/arpackSolver.pc.in "${PROJECT_BINARY_DIR}/EXAMPLES/MATRIX_MARKET/arpackSolver.pc" @ONLY)


install(TARGETS arpack
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(FILES "${PROJECT_BINARY_DIR}/SRC/arpack${LIBSUFFIX}.pc"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

if (MPI)
  install(TARGETS parpack
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  install(FILES "${PROJECT_BINARY_DIR}/PARPACK/SRC/MPI/parpack${LIBSUFFIX}.pc"
      DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
endif ()

if(ICB)
  install(FILES ICB/arpack.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/arpack)
  install(FILES ICB/arpack.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/arpack)
  if (MPI)
    install(FILES ICB/parpack.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/arpack)
    install(FILES ICB/parpack.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/arpack)
  endif()
  if (ICBEXMM)
    install(FILES EXAMPLES/MATRIX_MARKET/arpackSolver.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/arpack)
    install(FILES "${PROJECT_BINARY_DIR}/EXAMPLES/MATRIX_MARKET/arpackSolver.pc" DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
  endif()
endif()

install(FILES debug.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/arpack)
if(ICB)
  install(FILES ICB/debug_c.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/arpack)
  install(FILES ICB/debug_c.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/arpack)
endif()

install(FILES stat.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/arpack)
if(ICB)
  install(FILES ICB/stat_c.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/arpack)
  install(FILES ICB/stat_c.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/arpack)
endif()

configure_file(arpackdef.h.in "${PROJECT_BINARY_DIR}/arpackdef.h" @ONLY)
install(FILES "${PROJECT_BINARY_DIR}/arpackdef.h" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/arpack)

configure_file(arpackicb.h.in "${PROJECT_BINARY_DIR}/arpackicb.h" @ONLY)
install(FILES "${PROJECT_BINARY_DIR}/arpackicb.h" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/arpack)

# Provide find_package for arpack-ng to users.
configure_file(arpack-ng-config.cmake.in "${PROJECT_BINARY_DIR}/arpack-ng-config.cmake" @ONLY)
install(FILES "${PROJECT_BINARY_DIR}/arpack-ng-config.cmake" DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/arpack-ng) # find_package(arpack-ng)
configure_file(arpack-ng-config-version.cmake.in "${PROJECT_BINARY_DIR}/arpack-ng-config-version.cmake" @ONLY)
install(FILES "${PROJECT_BINARY_DIR}/arpack-ng-config-version.cmake" DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/arpack-ng)

# Packaging: ease arpack-ng distribution (precompiled binaries and sources tarballs).
set(CPACK_VERSION_MAJOR "${arpack_ng_MAJOR_VERSION}")
set(CPACK_VERSION_MINOR "${arpack_ng_MINOR_VERSION}")
set(CPACK_VERSION_PATCH "${arpack_ng_PATCH_VERSION}")
set(CPACK_SOURCE_IGNORE_FILES "/BUILD/" "/Build/" "/build/" "/local/") # Do not ship (eventual) build or install directories in tarballs.
set(CPACK_SOURCE_PACKAGE_FILE_NAME "arpack-ng-${CPACK_VERSION_MAJOR}.${CPACK_VERSION_MINOR}.${CPACK_VERSION_PATCH}")
include(CPack)

if (COVERALLS)
    set(arpack_TST_SRC
        ${arpack_SOURCE_DIR}/TESTS/bug_1323.f
        ${arpack_SOURCE_DIR}/TESTS/bug_142.f
        ${arpack_SOURCE_DIR}/TESTS/bug_142_gen.f
        ${arpack_SOURCE_DIR}/TESTS/bug_58_double.f
        ${arpack_SOURCE_DIR}/TESTS/bug_79_double_complex.f
        ${arpack_SOURCE_DIR}/TESTS/dnsimp.f
        ${arpack_SOURCE_DIR}/TESTS/mmio.f
        ${arpack_SOURCE_DIR}/TESTS/bug_1315_single.c
        ${arpack_SOURCE_DIR}/TESTS/bug_1315_double.c
        ${arpack_SOURCE_DIR}/TESTS/icb_arpack_c.c
        ${arpack_SOURCE_DIR}/TESTS/icb_arpack_cpp.cpp
       )
    file(GLOB_RECURSE arpack_EX_F_SRC "${arpack_SOURCE_DIR}/EXAMPLES/*/*.f")
    set(arpack_EX_CPP_SRC ${arpack_SOURCE_DIR}/EXAMPLES/MATRIX_MARKET/arpackmm.cpp)

    set(parpack_TST_SRC
        ${arpack_SOURCE_DIR}/PARPACK/TESTS/MPI/issue46.f
        ${arpack_SOURCE_DIR}/PARPACK/TESTS/MPI/icb_parpack_c.c
        ${arpack_SOURCE_DIR}/PARPACK/TESTS/MPI/icb_parpack_cpp.cpp
       )
    file(GLOB_RECURSE parpack_EX_F_SRC "${arpack_SOURCE_DIR}/PARPACK/EXAMPLES/MPI/*.f")

    # Create the coveralls target.
    # Also lists the c/cpp files for test purposes
    coveralls_setup(
        "${arpacksrc_STAT_SRCS} ${arpackutil_STAT_SRCS} ${arpacksrc_ICB} ${arpack_TST_SRC} ${arpack_EX_F_SRC} ${arpack_EX_CPP_SRC} ${parpacksrc_STAT_SRCS} ${parpackutil_STAT_SRCS} ${parpacksrc_ICB} ${parpack_TST_SRC} ${parpack_EX_F_SRC}" # The source files.
        ON                 # If we should upload.
        "${PROJECT_SOURCE_DIR}/cmake/") # (Optional) Alternate project cmake module path.
endif()

function(libsummary title include libraries)
    message("   -- ${title}:")
    foreach(inc ${include})
        message("      -- compile: ${inc}")
    endforeach()
    foreach(lib ${libraries})
        message("      -- link:    ${lib}")
    endforeach()
endfunction(libsummary)

function(cprsummary title compiler debug_flags minsizerel_flags release_flags relwithdebinfo_flags more_flags)
    message("   -- ${title}:      ${compiler}")
    if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug" OR "${CMAKE_BUILD_TYPE}" STREQUAL "DEBUG")
        message("   -- ${title}FLAGS: ${debug_flags} ${more_flags}")
    endif()
    if("${CMAKE_BUILD_TYPE}" STREQUAL "MinSizeRel" OR "${CMAKE_BUILD_TYPE}" STREQUAL "MINSIZEREL")
        message("   -- ${title}FLAGS: ${minsizerel_flags} ${more_flags}")
    endif()
    if("${CMAKE_BUILD_TYPE}" STREQUAL "Release" OR "${CMAKE_BUILD_TYPE}" STREQUAL "RELEASE")
        message("   -- ${title}FLAGS: ${release_flags} ${more_flags}")
    endif()
    if("${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo" OR "${CMAKE_BUILD_TYPE}" STREQUAL "RELWITHDEBINFO")
        message("   -- ${title}FLAGS: ${relwithdebinfo_flags} ${more_flags}")
    endif()
endfunction(cprsummary)

message("-- Configuration summary for arpack-ng-${arpack_ng_VERSION}:")
message("   -- prefix: ${CMAKE_INSTALL_PREFIX}")
message("   -- MPI: ${MPI}")
message("   -- ICB: ${ICB}")
message("   -- INTERFACE64: ${INTERFACE64}")
cprsummary("FC" "${CMAKE_Fortran_COMPILER}"
                "${CMAKE_Fortran_FLAGS_DEBUG}"
                "${CMAKE_Fortran_FLAGS_MINSIZEREL}"
                "${CMAKE_Fortran_FLAGS_RELEASE}"
                "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}"
                "${CMAKE_Fortran_FLAGS}")
if (ICB)
    cprsummary("CC" "${CMAKE_C_COMPILER}"
                    "${CMAKE_C_FLAGS_DEBUG}"
                    "${CMAKE_C_FLAGS_MINSIZEREL}"
                    "${CMAKE_C_FLAGS_RELEASE}"
                    "${CMAKE_C_FLAGS_RELWITHDEBINFO}"
                    "${CMAKE_C_FLAGS}")
    cprsummary("CXX" "${CMAKE_CXX_COMPILER}"
                     "${CMAKE_CXX_FLAGS_DEBUG}"
                     "${CMAKE_CXX_FLAGS_MINSIZEREL}"
                     "${CMAKE_CXX_FLAGS_RELEASE}"
                     "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}"
                     "${CMAKE_CXX_FLAGS}")
endif()
if (MPI)
    libsummary("MPIFC" "${MPI_Fortran_INCLUDE_DIRS}" "${MPI_Fortran_LIBRARIES}")
    if (ICB)
        libsummary("MPICC" "${MPI_C_INCLUDE_DIRS}" "${MPI_C_LIBRARIES}")
        libsummary("MPICXX" "${MPI_CXX_INCLUDE_DIRS}" "${MPI_CXX_LIBRARIES}")
    endif()
endif()
libsummary("BLAS" "" "${BLAS_LIBRARIES}")
libsummary("LAPACK" "" "${LAPACK_LIBRARIES}")
if (ICBEXMM)
    libsummary("EIGEN3" "${EIGEN3_INCLUDE_DIR}" "")
endif()
if (PYTHON3)
  libsummary("PYTHON" "${PYTHON_INCLUDE_DIRS}" "${PYTHON_LIBRARIES}")
  message("      -- exe:     ${PYTHON_EXECUTABLE}")
  libsummary("BOOST" "${Boost_INCLUDE_DIRS}" "${Boost_LIBRARIES}")
endif()
