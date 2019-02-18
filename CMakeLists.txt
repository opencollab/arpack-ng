cmake_minimum_required(VERSION 2.8.6)

if (NOT DEFINED CMAKE_BUILD_TYPE)
   set (CMAKE_BUILD_TYPE Release CACHE STRING "Build type")
endif ()

project(arpack C Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

# Adopted from https://github.com/feymark/arpack.git

if (POLICY CMP0042)
  # enable MACOSX_RPATH by default
  cmake_policy (SET CMP0042 NEW)
endif ()

option(MPI "Enable parallel support" OFF)
option(ICB "Enable support for *[ae]upd_c with ISO_C_BINDING" OFF)
option(ICBEXMM "Enable support for matrix market example based on ICB" OFF)
#option(SYSTEM_BLAS "Use system BLAS" ON)
#option(SYSTEM_LAPACK "Use system LAPACK" ON)
option(EXAMPLES "Compile ARPACK examples" OFF)
set(LIBSUFFIX ""
    CACHE STRING "suffix to add to ARPACK libraries names")
set(SYMBOLSUFFIX ""
    CACHE STRING "suffix to add to ARPACK, BLAS and LAPACK function names")
option(INTERFACE64 "use the 64-bit integer interface (ILP64) for ARPACK, BLAS and LAPACK")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

if (COVERALLS)
    include(Coveralls)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -O0 -fprofile-arcs -ftest-coverage")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O0 -fprofile-arcs -ftest-coverage")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -fprofile-arcs -ftest-coverage")
    # The no space is by design: issue in cmake. See CMP0004.
    set(EXTRA_LDLAGS "${EXTRA_LDLAGS}-lgcov")
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
        target_link_libraries(${lwe} arpack ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${EXTRA_LDLAGS})
        add_test(NAME "${lwe}_ex" COMMAND ${lwe} WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
    endforeach()
endfunction(examples)

function(pexamples list_name)
    foreach(l ${${list_name}})
        get_filename_component(lwe ${l} NAME_WE)
        add_executable(${lwe} ${parpackexample_DIR}/${l} )
        target_link_libraries(${lwe} parpack arpack  ${MPI_Fortran_LIBRARIES})
        add_test(NAME "${lwe}_ex" COMMAND mpirun -n 2 ${lwe} WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
    endforeach()
endfunction(pexamples)

if (ICBEXMM)
    find_package(Eigen3)
    if (NOT Eigen3_FOUND) # If not found, piggy-back pkg-config files.
        message(WARNING "CMake didn't find the Eigen3 package. Try to look for pkg-config file...")
        find_package(PkgConfig REQUIRED)
        pkg_check_modules(EIGEN3 REQUIRED eigen3>=3.2)
    endif ()
    find_program (BASH_PROGRAM bash)

    # Look for headers.
    find_path(EIGEN3_SPARSE_DIR NAMES Sparse PATHS ${EIGEN3_INCLUDE_DIRS} PATH_SUFFIXES Eigen)
    if (NOT EIGEN3_SPARSE_DIR)
        message(FATAL_ERROR "-- Eigen/Sparse header not found.")
    endif ()
    find_path(EIGEN3_ITERATIVE_SOLVER_DIR NAMES IterativeLinearSolvers PATHS ${EIGEN3_INCLUDE_DIRS} PATH_SUFFIXES Eigen)
    if (NOT EIGEN3_ITERATIVE_SOLVER_DIR)
        message(FATAL_ERROR "-- Eigen/IterativeLinearSolvers header not found.")
    endif ()
    find_path(EIGEN3_SPARSELU_SOLVER_DIR NAMES SparseLU PATHS ${EIGEN3_INCLUDE_DIRS} PATH_SUFFIXES Eigen)
    if (NOT EIGEN3_SPARSELU_SOLVER_DIR)
        message(FATAL_ERROR "-- Eigen/SparseLU header not found.")
    endif ()
    find_path(EIGEN3_SPARSEQR_SOLVER_DIR NAMES SparseQR PATHS ${EIGEN3_INCLUDE_DIRS} PATH_SUFFIXES Eigen)
    if (NOT EIGEN3_SPARSEQR_SOLVER_DIR)
        message(FATAL_ERROR "-- Eigen/SparseQR header not found.")
    endif ()
    find_path(EIGEN3_SPARSECHOLESKY_SOLVER_DIR NAMES SparseCholesky PATHS ${EIGEN3_INCLUDE_DIRS} PATH_SUFFIXES Eigen)
    if (NOT EIGEN3_SPARSECHOLESKY_SOLVER_DIR)
        message(FATAL_ERROR "-- Eigen/SparseCholesky header not found.")
    endif ()
    find_path(EIGEN3_DENSE_DIR NAMES Dense PATHS ${EIGEN3_INCLUDE_DIRS} PATH_SUFFIXES Eigen)
    if (NOT EIGEN3_DENSE_DIR)
        message(FATAL_ERROR "-- Eigen/Dense header not found.")
    endif ()
    find_path(EIGEN3_DENSELU_DIR NAMES LU PATHS ${EIGEN3_INCLUDE_DIRS} PATH_SUFFIXES Eigen)
    if (NOT EIGEN3_DENSELU_DIR)
        message(FATAL_ERROR "-- Eigen/LU header not found.")
    endif ()
    find_path(EIGEN3_DENSEQR_DIR NAMES QR PATHS ${EIGEN3_INCLUDE_DIRS} PATH_SUFFIXES Eigen)
    if (NOT EIGEN3_DENSEQR_DIR)
        message(FATAL_ERROR "-- Eigen/QR header not found.")
    endif ()
    find_path(EIGEN3_DENSECHOLESKY_DIR NAMES Cholesky PATHS ${EIGEN3_INCLUDE_DIRS} PATH_SUFFIXES Eigen)
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
    # ICB saves you from old-fashion-boring-cumbersome-fortran/C crap...

    if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -ffixed-line-length-none")
    elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpp -extend-source")
    else ()
      message(WARNING "build script does not know how to enable your Fortran compiler's preprocessor and support for lines longer than 72 characters: set them manually via FFLAGS.")
    endif ()

    set(SCALARFUNS naitr napps naup2 naupd nconv neigh neupd ngets statn saitr sapps saup2 saupd sconv seigt seupd sgets stats getv0 sortc sortr sesrt stqrb)
    set(COMPLEXFUNS, naitr napps naup2 naupd neigh neupd ngets statn getv0 sortc)

    set(BLASFUNS1 axpy copy gemv geqr2 lacpy lae2 lahqr lanhs larnv lartg lascl laset lasrt scal trevc trmm trsen gbmv gbtrf gbtrs gttrf gttrs pttrf pttrs)
    set(BLASFUNS2 dot ger labad laev2 lamch lanst lanv2 lapy2 larf larfg lasr nrm2 orm2r rot steqr swap)
    set(BLASFUNS3 dotc geru unm2r)
    set(BLASFUNS4 COPY LABAD LAMCH LANHS LANV2 LARFG ROT GEMV)
    set(BLASFUNS5 scnrm2 dznrm2 csscal zdscal)

    if (SYMBOLSUFFIX)
        foreach (f IN LISTS SCALARFUNS BLASFUNS1 BLASFUNS2)
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Ds${f}=s${f}${SYMBOLSUFFIX} -Dd${f}=d${f}${SYMBOLSUFFIX}")
        endforeach ()

        foreach (f IN LISTS COMPLEXFUNS BLASFUNS1 BLASFUNS3)
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Dc${f}=c${f}${SYMBOLSUFFIX} -Dz${f}=z${f}${SYMBOLSUFFIX}")
        endforeach ()

        foreach (f IN LISTS BLASFUNS4)
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DS${f}=S${f}${SYMBOLSUFFIX} -DD${f}=D${f}${SYMBOLSUFFIX}")
        endforeach ()

        foreach (f IN LISTS BLASFUNS5)
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -D${f}=${f}${SYMBOLSUFFIX}")
        endforeach ()
    endif ()

    set(CFUNS sgemm snaupd sneupd dnaupd dneupd cheev)
    foreach (f IN LISTS CFUNS)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D${f}=${f}${SYMBOLSUFFIX}")
        list(APPEND CFUNS_SUFFIXED ${f}${SYMBOLSUFFIX})
    endforeach ()

    include(FortranCInterface)
    FortranCInterface_HEADER(FCMangle.h SYMBOLS ${CFUNS_SUFFIXED})

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DINCLUDE_FCMANGLE")

    FortranCInterface_VERIFY()
endif ()

#if (SYSTEM_BLAS)
    find_package(BLAS REQUIRED)
#endif()
#if (BLAS_LIBRARIES)
#    set(SYSTEM_BLAS ON)
#endif()

if (MPI)
    include(FindMPI)
    find_package(MPI REQUIRED)
    include_directories(${MPI_Fortran_INCLUDE_PATH})
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAG}")

    if(ICB)
        include_directories(${MPI_C_INCLUDE_PATH})
        include_directories(${MPI_CXX_INCLUDE_PATH})
        include(CheckSymbolExists)
        check_symbol_exists(MPI_Comm_c2f "${MPI_C_INCLUDE_PATH}/mpi.h" MPI_Comm_c2f_FOUND)
        if(NOT ${MPI_Comm_c2f_FOUND})
          message(FATAL_ERROR "symbol MPI_Comm_c2f does not exist")
        endif()
    endif()
endif()

#if (SYSTEM_LAPACK)
    find_package(LAPACK REQUIRED)
#endif()
#if (BLAS_LIBRARIES)
#    set(SYSTEM_LAPACK ON)
#endif()

if (MPI)
    set(parpackutil_DIR ${arpack_SOURCE_DIR}/PARPACK/UTIL/)
    set(parpacksrc_DIR ${arpack_SOURCE_DIR}/PARPACK/SRC/)
endif()

#if (NOT SYSTEM_BLAS)
#    file(GLOB arpackblas_STAT_SRCS blas/*.f)
#endif()
#if (NOT SYSTEM_LAPACK)
#    file(GLOB arpacklapack_STAT_SRCS lapack/*.f)
#endif()
#file(GLOB arpackutil_STAT_SRCS util/*.f)
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

if(ICB)
  install(FILES ICB/arpack.h DESTINATION include/arpack)
  install(FILES ICB/arpack.hpp DESTINATION include/arpack)
  if (MPI)
    install(FILES ICB/parpack.h DESTINATION include/arpack)
    install(FILES ICB/parpack.hpp DESTINATION include/arpack)
  endif()
endif()

install(FILES debug.h DESTINATION include/arpack)
if(ICB)
  install(FILES ICB/debug_c.h DESTINATION include/arpack)
  install(FILES ICB/debug_c.hpp DESTINATION include/arpack)
endif()

install(FILES stat.h DESTINATION include/arpack)
if(ICB)
  install(FILES ICB/stat_c.h DESTINATION include/arpack)
  install(FILES ICB/stat_c.hpp DESTINATION include/arpack)
endif()

target_link_libraries(arpack ${BLAS_LIBRARIES})
target_link_libraries(arpack ${LAPACK_LIBRARIES} ${EXTRA_LDLAGS})
set_target_properties(arpack PROPERTIES OUTPUT_NAME arpack${LIBSUFFIX})
set_target_properties(arpack PROPERTIES VERSION 2.1.0)
set_target_properties(arpack PROPERTIES SOVERSION 2)

if (MPI)
    # use -DBUILD_SHARED_LIBS=ON|OFF to control static/shared
    add_library(parpack ${parpacksrc_STAT_SRCS} ${parpackutil_STAT_SRCS} ${parpacksrc_ICB})

    if (ICB)
        target_include_directories(parpack PUBLIC ${MPI_C_INCLUDE_DIRS})
        target_include_directories(parpack PUBLIC ${MPI_CXX_INCLUDE_DIRS})
        target_link_libraries(parpack ${MPI_C_LIBRARIES})
        target_link_libraries(parpack ${MPI_CXX_LIBRARIES})
    endif()
    target_include_directories(parpack PUBLIC ${MPI_Fortran_INCLUDE_DIRS})
    target_link_libraries(parpack ${MPI_Fortran_LIBRARIES})
    target_link_libraries(parpack arpack)
    set_target_properties(parpack PROPERTIES OUTPUT_NAME parpack${LIBSUFFIX})
    set_target_properties(parpack PROPERTIES VERSION 2.1.0)
    set_target_properties(parpack PROPERTIES SOVERSION 2)
endif ()

if (INTERFACE64)
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

message("-- Generating arpackdef.h")
configure_file(arpackdef.h.in "${PROJECT_BINARY_DIR}/arpackdef.h" @ONLY)
include_directories("${PROJECT_BINARY_DIR}") # Find arpackdef.h
install(FILES "${PROJECT_BINARY_DIR}/arpackdef.h" DESTINATION include/arpack)

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
target_link_libraries(dnsimp_test arpack ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${EXTRA_LDLAGS})
add_custom_command(TARGET dnsimp_test POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/TESTS/testA.mtx testA.mtx
)
add_test(dnsimp_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/dnsimp)

if (NOT ICB)
    add_executable(bug_1315_single TESTS/bug_1315_single.c)
    target_link_libraries(bug_1315_single arpack ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${EXTRA_LDLAGS})
    add_test(bug_1315_single_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/bug_1315_single)

    add_executable(bug_1315_double TESTS/bug_1315_double.c)
    target_link_libraries(bug_1315_double arpack ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${EXTRA_LDLAGS})
    add_test(bug_1315_double_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/bug_1315_double)
endif()

add_executable(bug_1323 TESTS/bug_1323.f)
target_link_libraries(bug_1323 arpack ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${EXTRA_LDLAGS})
add_test(bug_1323_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/bug_1323)

add_executable(bug_58_double TESTS/bug_58_double.f)
target_link_libraries(bug_58_double arpack ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${EXTRA_LDLAGS})
add_test(bug_58_double_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/bug_58_double)

add_executable(bug_79_double_complex TESTS/bug_79_double_complex.f)
target_link_libraries(bug_79_double_complex arpack ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${EXTRA_LDLAGS})
add_test(bug_79_double_complex_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/bug_79_double_complex)

add_executable(bug_142 TESTS/bug_142.f)
target_link_libraries(bug_142 arpack ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${EXTRA_LDLAGS})
add_test(bug_142_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/bug_142)

add_executable(bug_142_gen TESTS/bug_142_gen.f)
target_link_libraries(bug_142_gen arpack ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${EXTRA_LDLAGS})
add_test(bug_142_gen_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/bug_142_gen)

if(MPI)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/PARPACK/TESTS/MPI)

  add_executable(issue46 PARPACK/TESTS/MPI/issue46.f)
  target_link_libraries(issue46 parpack arpack ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${EXTRA_LDLAGS})
  add_test(issue46_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/issue46)
endif()

if(ICB)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/TESTS)

  add_executable(icb_arpack_c TESTS/icb_arpack_c.c)
  target_include_directories(icb_arpack_c PUBLIC ${PROJECT_SOURCE_DIR}/ICB) # Get arpack.h
  target_link_libraries(icb_arpack_c arpack ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${EXTRA_LDLAGS})
  add_test(icb_arpack_c_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/icb_arpack_c)

  add_executable(icb_arpack_cpp TESTS/icb_arpack_cpp.cpp)
  target_include_directories(icb_arpack_cpp PUBLIC ${PROJECT_SOURCE_DIR}/ICB) # Get arpack.hpp
  target_link_libraries(icb_arpack_cpp arpack ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${EXTRA_LDLAGS})
  add_test(icb_arpack_cpp_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/icb_arpack_cpp)

  if (ICBEXMM)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/EXAMPLES/MATRIX_MARKET)

    add_executable(arpackmm EXAMPLES/MATRIX_MARKET/arpackmm.cpp)
    target_include_directories(arpackmm PUBLIC ${PROJECT_SOURCE_DIR}/ICB ${EIGEN3_INCLUDE_DIRS}) # Get arpack.h + eigen
    target_link_libraries(arpackmm arpack ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${EXTRA_LDLAGS})
    configure_file(EXAMPLES/MATRIX_MARKET/As.mtx ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/As.mtx)
    configure_file(EXAMPLES/MATRIX_MARKET/An.mtx ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/An.mtx)
    configure_file(EXAMPLES/MATRIX_MARKET/Az.mtx ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Az.mtx)
    configure_file(EXAMPLES/MATRIX_MARKET/B.mtx  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/B.mtx)
    configure_file(EXAMPLES/MATRIX_MARKET/Bz.mtx ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Bz.mtx)
    configure_file(EXAMPLES/MATRIX_MARKET/arpackmm.sh ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/arpackmm.sh)
    add_test(NAME arpackmm_tst WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} COMMAND ${BASH_PROGRAM} arpackmm.sh)
  endif()

  if (MPI)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/PARPACK/TESTS/MPI)

    add_executable(icb_parpack_c PARPACK/TESTS/MPI/icb_parpack_c.c)
    target_include_directories(icb_parpack_c PUBLIC ${PROJECT_SOURCE_DIR}/ICB ${MPI_C_INCLUDE_DIRS}) # Get parpack.h mpi.h
    target_link_libraries(icb_parpack_c parpack arpack ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${EXTRA_LDLAGS} ${MPI_C_LIBRARIES})
    add_test(icb_parpack_c_tst mpirun -n 2 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/icb_parpack_c)

    add_executable(icb_parpack_cpp PARPACK/TESTS/MPI/icb_parpack_cpp.cpp)
    target_include_directories(icb_parpack_cpp PUBLIC ${PROJECT_SOURCE_DIR}/ICB ${MPI_CXX_INCLUDE_DIRS}) # Get parpack.hpp mpi.h
    target_link_libraries(icb_parpack_cpp parpack arpack ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${EXTRA_LDLAGS} ${MPI_CXX_LIBRARIES})
    add_test(icb_parpack_cpp_tst mpirun -n 2 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/icb_parpack_cpp)
  endif()
endif()

############################
# install
############################
# 'make install' to the correct location
include(GNUInstallDirs)

install(TARGETS arpack
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

if (MPI)
  install(TARGETS parpack
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif ()

# Provide find_package for arpack-ng to users.
configure_file(arpack-ng-config.cmake.in "${PROJECT_BINARY_DIR}/arpack-ng-config.cmake" @ONLY)
install(FILES "${PROJECT_BINARY_DIR}/arpack-ng-config.cmake" DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake") # find_package(arpack-ng)
set(arpack_ng_MAJOR_VERSION 3)
set(arpack_ng_MINOR_VERSION 7)
set(arpack_ng_PATCH_VERSION 0)
set(arpack_ng_VERSION ${arpack_ng_MAJOR_VERSION}.${arpack_ng_MINOR_VERSION}.${arpack_ng_PATCH_VERSION})
configure_file(arpack-ng-config-version.cmake.in "${PROJECT_BINARY_DIR}/arpack-ng-config-version.cmake" @ONLY)
install(FILES "${PROJECT_BINARY_DIR}/arpack-ng-config-version.cmake" DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake")

# Packaging: ease arpack-ng distribution (precompiled binaries and sources tarballs).
set(CPACK_VERSION_MAJOR "${arpack_ng_MAJOR_VERSION}")
set(CPACK_VERSION_MINOR "${arpack_ng_MINOR_VERSION}")
set(CPACK_VERSION_PATCH "${arpack_ng_PATCH_VERSION}")
set(CPACK_SOURCE_IGNORE_FILES "/BUILD/" "/Build/" "/build/" "/local/") # Do not ship (eventual) build or install directories in tarballs.
set(CPACK_SOURCE_PACKAGE_FILE_NAME "arpack-ng-${CPACK_VERSION_MAJOR}.${CPACK_VERSION_MINOR}.${CPACK_VERSION_PATCH}")
include(CPack)

# Convert variable names to those expected by the .pc file.
set(prefix ${CMAKE_INSTALL_PREFIX})
set(exec_prefix \${prefix})
set(libdir \${exec_prefix}/${CMAKE_INSTALL_LIBDIR})
set(includedir \${prefix}/${CMAKE_INSTALL_INCLUDEDIR})
set(PACKAGE_NAME ${PROJECT_NAME})
set(PACKAGE_VERSION ${arpack_ng_VERSION})
set(PACKAGE_URL "https://github.com/opencollab/arpack-ng/")
set(LAPACK_LIBS ${LAPACK_LIBRARIES})
set(BLAS_LIBS ${BLAS_LIBRARIES})
configure_file(arpack.pc.in arpack${LIBSUFFIX}.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/arpack${LIBSUFFIX}.pc
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

# We don't want this to run on every build.
option(COVERALLS "Generate coveralls data" OFF)


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
    libsummary("MPIFC" "${MPI_Fortran_INCLUDE_PATH}" "${MPI_Fortran_LIBRARIES}")
    if (ICB)
        libsummary("MPICC" "${MPI_C_INCLUDE_PATH}" "${MPI_C_LIBRARIES}")
        libsummary("MPICXX" "${MPI_CXX_INCLUDE_PATH}" "${MPI_CXX_LIBRARIES}")
    endif()
endif()
libsummary("BLAS" "" "${BLAS_LIBRARIES}")
libsummary("LAPACK" "" "${LAPACK_LIBRARIES}")
if (ICBEXMM)
    libsummary("EIGEN3" "${EIGEN3_INCLUDE_DIRS}" "")
endif()
