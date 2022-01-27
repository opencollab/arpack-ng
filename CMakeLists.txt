cmake_minimum_required(VERSION 3.0)

if (NOT DEFINED CMAKE_BUILD_TYPE)
   set (CMAKE_BUILD_TYPE Release CACHE STRING "Build type")
endif ()

project(arpack C Fortran)

set(arpack_target ARPACK)
set(parpack_target PARPACK)
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

# Suffixes: LIBSUFFIX modify ONLY libraries names, ITF64SUFFIX modify BOTH libraries AND include directory names.
set(LIBSUFFIX ""
    CACHE STRING "suffix to add to ARPACK libraries names")
set(ITF64SUFFIX ""
    CACHE STRING "suffix to add to ARPACK include directory and libraries names (use with INTERFACE64)")
set(SYMBOLSUFFIX ""
    CACHE STRING "suffix to add to ARPACK, BLAS and LAPACK function names")
option(INTERFACE64 "use the 64-bit integer interface (ILP64) for ARPACK, BLAS and LAPACK")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Utility functions for creating tests and detecting features
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include(ArpackFunctions)

# 'make install' to the correct location
include(GNUInstallDirs)
set(ARPACK_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}/arpack-ng${ITF64SUFFIX}")
set(ARPACK_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/arpackng${LIBSUFFIX}${ITF64SUFFIX}")

if (PYTHON3)
    enable_language(C CXX) # Boost requirement.
    set(CMAKE_CXX_STANDARD 14) # Boost requirement.

    find_package(PythonInterp 3 REQUIRED)
    find_package(PythonLibs 3 REQUIRED)
    find_package(Boost 1.78 COMPONENTS python${BOOST_PYTHON_LIBSUFFIX} numpy${BOOST_PYTHON_LIBSUFFIX} REQUIRED)

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

    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAG}")

    # Check if we can use ISO_C_BINDING provided by MPI.
    file(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/PROG_ICB.f90
         "
          PROGRAM PROG_ICB
          USE :: mpi_f08
          IMPLICIT NONE
          type(MPI_Comm)    comm
          type(MPI_Status)  status
          END PROGRAM PROG_ICB
         ")
         try_compile(COMPILE_ICB ${CMAKE_BINARY_DIR} ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/PROG_ICB.f90 LINK_LIBRARIES MPI::MPI_Fortran)
    if(NOT ${COMPILE_ICB})
      message("-- MPI library does not support iso_c_binding.")
      set(HAVE_MPI_ICB 0)
    else()
      message("-- MPI library does support iso_c_binding.")
      set(HAVE_MPI_ICB 1)
      add_compile_definitions(HAVE_MPI_ICB=1)
    endif()

    # As MPI can be used with or without ISO_C_BINDING (#ifdef), we need to preprocess code before compiling.
    if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")
    elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpp")
    else ()
      message(WARNING "build script does not know how to preprocess Fortran code: set it manually via FFLAGS.")
    endif ()

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

        if (NOT TARGET MPI::MPI_CXX) # Search only if not already found by upper CMakeLists.txt
            include(FindMPI)
            find_package(MPI REQUIRED COMPONENTS CXX)

            if (NOT TARGET MPI::MPI_CXX) # Create target "at hand" to ensure compatibility if cmake version < 3.9
                add_library(MPI::MPI_CXX INTERFACE IMPORTED)
                set_target_properties(MPI::MPI_CXX PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${MPI_CXX_INCLUDE_DIRS}")
                set_target_properties(MPI::MPI_CXX PROPERTIES INTERFACE_LINK_LIBRARIES      "${MPI_CXX_LIBRARIES}")
            endif()
        endif()

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

# Find LAPACK

if (NOT TARGET LAPACK::LAPACK) # Search only if not already found by upper CMakeLists.txt
    find_package(LAPACK REQUIRED)

    # LAPACK::LAPACK target was already created at this point by FindLAPACK.cmake if cmake version >= 3.18
    if (NOT TARGET LAPACK::LAPACK) # Create target "at hand" to ensure compatibility if cmake version < 3.18
        add_library(LAPACK::LAPACK INTERFACE IMPORTED)
        set_target_properties(LAPACK::LAPACK PROPERTIES INTERFACE_LINK_LIBRARIES "${LAPACK_LIBRARIES}")
    endif()
endif()

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
add_library(${arpack_target} ${arpackutil_STAT_SRCS} ${arpacksrc_STAT_SRCS} ${arpacksrc_ICB})

target_link_libraries(${arpack_target}
  PUBLIC
  $<INSTALL_INTERFACE:$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,,BLAS::BLAS>>
  $<INSTALL_INTERFACE:$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,,LAPACK::LAPACK>>
  $<BUILD_INTERFACE:BLAS::BLAS>
  $<BUILD_INTERFACE:LAPACK::LAPACK>
)
target_link_options(${arpack_target}
  PUBLIC
  "$<INSTALL_INTERFACE:${EXTRA_LDFLAGS}>"
  "$<BUILD_INTERFACE:${EXTRA_LDFLAGS}>"
)
set_target_properties(${arpack_target} PROPERTIES OUTPUT_NAME arpack${LIBSUFFIX}${ITF64SUFFIX})
set_target_properties(${arpack_target} PROPERTIES VERSION 2.1.0)
set_target_properties(${arpack_target} PROPERTIES SOVERSION 2)

target_include_directories(${arpack_target}
  PUBLIC
  # Exported location of headers
  $<INSTALL_INTERFACE:${ARPACK_INSTALL_INCLUDEDIR}>
  # Find arpackdef.h, arpackicb.h, stat*.h, debug*.h at build time
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  # For ICB interface
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/ICB>
)
if (ICBEXMM)
	target_include_directories(${arpack_target} PUBLIC ${EIGEN3_INCLUDE_DIR})
endif()

if (MPI)
    # use -DBUILD_SHARED_LIBS=ON|OFF to control static/shared
    add_library(${parpack_target} ${parpacksrc_STAT_SRCS} ${parpackutil_STAT_SRCS} ${parpacksrc_ICB})
    set_target_properties(${parpack_target} PROPERTIES OUTPUT_NAME parpack${LIBSUFFIX}${ITF64SUFFIX})
    set_target_properties(${parpack_target} PROPERTIES VERSION 2.1.0)
    set_target_properties(${parpack_target} PROPERTIES SOVERSION 2)
	target_link_libraries(${parpack_target}
	  PUBLIC
	  ${arpack_target}
	  $<INSTALL_INTERFACE:$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,,MPI::MPI_Fortran>>
	  $<BUILD_INTERFACE:MPI::MPI_Fortran>
	)
endif()

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

    # TODO: this needs full re-write of parpack to support ILP64...
    if (MPI)
      message(FATAL_ERROR "Parallel arpack does not support ILP64.")
    endif()
else ()
    set(INTERFACE64 0)
endif ()

if (PYTHON3)
    python_add_module(pyarpack ${arpackutil_STAT_SRCS} ${arpacksrc_STAT_SRCS} ${arpacksrc_ICB} ${PROJECT_SOURCE_DIR}/EXAMPLES/PYARPACK/pyarpack.cpp)
	target_include_directories(pyarpack
		PRIVATE
		${PROJECT_BINARY_DIR}
		${PROJECT_SOURCE_DIR}
		${PROJECT_SOURCE_DIR}/ICB
		${PROJECT_SOURCE_DIR}/EXAMPLES/MATRIX_MARKET
		${PROJECT_SOURCE_DIR}/EXAMPLES/PYARPACK
		${EIGEN3_INCLUDE_DIR}
		${Boost_INCLUDE_DIRS}
		${PYTHON_INCLUDE_DIRS})
	target_link_libraries(pyarpack
		BLAS::BLAS LAPACK::LAPACK ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
	install(TARGETS pyarpack
		ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/pyarpack
		LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/pyarpack)
endif()

############################
# EXAMPLES AND TESTING
############################

enable_testing()
set(CMAKE_CTEST_COMMAND ctest -V)

if (EXAMPLES)
	add_subdirectory("EXAMPLES")
	add_subdirectory("PARPACK/EXAMPLES/MPI")
endif()
add_subdirectory("TESTS")
add_subdirectory("PARPACK/TESTS/MPI")
add_subdirectory("EXAMPLES/PYARPACK")
add_subdirectory("EXAMPLES/MATRIX_MARKET")

############################
# install
############################

# Convert variable names to those expected by the .pc file.
set(prefix ${CMAKE_INSTALL_PREFIX})
set(exec_prefix \${prefix})
set(libdir ${CMAKE_INSTALL_FULL_LIBDIR})
set(includedir ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR})
set(PACKAGE_NAME ${PROJECT_NAME})
set(PACKAGE_VERSION ${arpack_ng_VERSION})
set(PACKAGE_URL "https://github.com/opencollab/arpack-ng/")

# Convert (LAPACK|BLAS)_LIBRARIES to -l flags, store in ARPACK_PC_LIBS_PRIVATE.
# LAPACK_LIBRARIES contains transitive deps, no need to parse BLAS_LIBRARIES.
set(ARPACK_PC_LIBS_PRIVATE)
foreach(lib ${LAPACK_LIBRARIES})
    get_filename_component(libname ${lib} NAME)
    string(REGEX REPLACE "^(lib)?([^.]+).*$" "-l\\2" libname ${libname})
    list(APPEND ARPACK_PC_LIBS_PRIVATE "${libname}")
endforeach()
string(REPLACE ";" " " ARPACK_PC_LIBS_PRIVATE "${ARPACK_PC_LIBS_PRIVATE}")

set(PARPACK_PC_LIBS_PRIVATE)
foreach(lib ${LAPACK_LIBRARIES} ${MPI_Fortran_LIBRARIES})
    get_filename_component(libname ${lib} NAME)
    string(REGEX REPLACE "^(lib)?([^.]+).*$" "-l\\2" libname ${libname})
    list(APPEND PARPACK_PC_LIBS_PRIVATE "${libname}")
endforeach()
string(REPLACE ";" " " PARPACK_PC_LIBS_PRIVATE "${PARPACK_PC_LIBS_PRIVATE}")

configure_file(pkg-config/arpack.pc.in "${PROJECT_BINARY_DIR}/arpack${LIBSUFFIX}${ITF64SUFFIX}.pc" @ONLY)
configure_file(pkg-config/parpack.pc.in "${PROJECT_BINARY_DIR}/parpack${LIBSUFFIX}${ITF64SUFFIX}.pc" @ONLY)
configure_file(pkg-config/arpackSolver.pc.in "${PROJECT_BINARY_DIR}/arpackSolver${LIBSUFFIX}${ITF64SUFFIX}.pc" @ONLY)


install(TARGETS ${arpack_target}
    EXPORT arpackngTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(FILES "${PROJECT_BINARY_DIR}/arpack${LIBSUFFIX}${ITF64SUFFIX}.pc"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

if (MPI)
  install(TARGETS ${parpack_target}
      EXPORT arpackngTargets
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  install(FILES "${PROJECT_BINARY_DIR}/parpack${LIBSUFFIX}${ITF64SUFFIX}.pc"
      DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
endif ()

if(ICB)
  install(FILES ICB/arpack.h DESTINATION "${ARPACK_INSTALL_INCLUDEDIR}")
  install(FILES ICB/arpack.hpp DESTINATION ${ARPACK_INSTALL_INCLUDEDIR})
  if (MPI)
    install(FILES ICB/parpack.h DESTINATION ${ARPACK_INSTALL_INCLUDEDIR})
    install(FILES ICB/parpack.hpp DESTINATION ${ARPACK_INSTALL_INCLUDEDIR})
  endif()
  if (ICBEXMM)
    install(FILES EXAMPLES/MATRIX_MARKET/arpackSolver.hpp DESTINATION ${ARPACK_INSTALL_INCLUDEDIR})
    install(FILES "${PROJECT_BINARY_DIR}/arpackSolver${LIBSUFFIX}${ITF64SUFFIX}.pc" DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
  endif()
endif()

install(FILES debug.h DESTINATION ${ARPACK_INSTALL_INCLUDEDIR})
install(FILES debugF90.h DESTINATION ${ARPACK_INSTALL_INCLUDEDIR})
if(ICB)
  install(FILES ICB/debug_c.h DESTINATION ${ARPACK_INSTALL_INCLUDEDIR})
  install(FILES ICB/debug_c.hpp DESTINATION ${ARPACK_INSTALL_INCLUDEDIR})
endif()

install(FILES stat.h DESTINATION ${ARPACK_INSTALL_INCLUDEDIR})
install(FILES statF90.h DESTINATION ${ARPACK_INSTALL_INCLUDEDIR})
if(ICB)
  install(FILES ICB/stat_c.h DESTINATION ${ARPACK_INSTALL_INCLUDEDIR})
  install(FILES ICB/stat_c.hpp DESTINATION ${ARPACK_INSTALL_INCLUDEDIR})
endif()

configure_file(arpackdef.h.in "${PROJECT_BINARY_DIR}/arpackdef.h" @ONLY)
install(FILES "${PROJECT_BINARY_DIR}/arpackdef.h" DESTINATION ${ARPACK_INSTALL_INCLUDEDIR})

configure_file(arpackicb.h.in "${PROJECT_BINARY_DIR}/arpackicb.h" @ONLY)
install(FILES "${PROJECT_BINARY_DIR}/arpackicb.h" DESTINATION ${ARPACK_INSTALL_INCLUDEDIR})

# Provide find_package for arpack-ng to users.
install(EXPORT arpackngTargets
  NAMESPACE ARPACK::
  DESTINATION "${ARPACK_INSTALL_CMAKEDIR}"
  )
configure_file(cmake/arpackng-config.cmake.in "${PROJECT_BINARY_DIR}/arpackng-config.cmake" @ONLY)
configure_file(cmake/arpackng-config-version.cmake.in "${PROJECT_BINARY_DIR}/arpackng-config-version.cmake" @ONLY)
install(FILES "${PROJECT_BINARY_DIR}/arpackng-config.cmake"
              "${PROJECT_BINARY_DIR}/arpackng-config-version.cmake"
  DESTINATION "${ARPACK_INSTALL_CMAKEDIR}") # find_package(arpackng)

configure_file(cmake/tstCMakeInstall.sh.in ${PROJECT_BINARY_DIR}/tstCMakeInstall.sh @ONLY)

# Packaging: ease arpack-ng distribution (precompiled binaries and sources tarballs).
set(CPACK_VERSION_MAJOR "${arpack_ng_MAJOR_VERSION}")
set(CPACK_VERSION_MINOR "${arpack_ng_MINOR_VERSION}")
set(CPACK_VERSION_PATCH "${arpack_ng_PATCH_VERSION}")
set(CPACK_SOURCE_IGNORE_FILES "/BUILD/" "/Build/" "/build/" "/local/") # Do not ship (eventual) build or install directories in tarballs.
set(CPACK_SOURCE_PACKAGE_FILE_NAME "arpack-ng-${CPACK_VERSION_MAJOR}.${CPACK_VERSION_MINOR}.${CPACK_VERSION_PATCH}")
include(CPack)

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
message("   -- MPI: ${MPI} (ICB provided ${HAVE_MPI_ICB})")
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
