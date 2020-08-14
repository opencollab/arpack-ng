# Install arpack-ng in CMAKE_BINARY_DIR.
execute_process(COMMAND make install DESTDIR=${CMAKE_BINARY_DIR}/local COMMAND_ECHO STDOUT RESULT_VARIABLE CMD_RESULT)
if (CMD_RESULT)
    message(FATAL_ERROR "Error in TestFindPackage - install")
endif ()

# Use find_package: call cmake to use this install of arpack-ng.
set(CMAKE_PREFIX_PATH "${CMAKE_BINARY_DIR}/local/usr/local/lib/cmake") # Add directory that contains arpack-ng-config.cmake to CMAKE_PREFIX_PATH.
find_package(arpack-ng CONFIG REQUIRED) # Use CONFIG as we look for arpack-ng-config.cmake (not FindXXX.cmake).

# Build with arpack.
add_executable(dnsimp_find_package ${PROJECT_SOURCE_DIR}/EXAMPLES/SIMPLE/dnsimp.f)
target_include_directories(dnsimp_find_package PUBLIC ARPACK::ARPACK)
target_link_libraries(dnsimp_find_package ARPACK::ARPACK)

# Build with parpack.
if (MPI)
    add_executable(pdndrv1_find_package ${PROJECT_SOURCE_DIR}/PARPACK/EXAMPLES/MPI/pdndrv1.f)
    target_include_directories(pdndrv1_find_package PUBLIC PARPACK::PARPACK)
    target_link_libraries(pdndrv1_find_package PARPACK::PARPACK)
endif ()
