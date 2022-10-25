# Config file for the arpack-ng package.
#
# To use arpack from CMake, use ARPACK::ARPACK target:
#   find_package(arpackng)
#   add_executable(main main.f)
#   target_include_directories(main INTERFACE ARPACK::ARPACK)
#   target_link_libraries(main ARPACK::ARPACK)
#
# To use parpack from CMake, use PARPACK::PARPACK target:
#   find_package(arpackng)
#   add_executable(main main.f)
#   target_include_directories(main INTERFACE PARPACK::PARPACK)
#   target_link_libraries(main PARPACK::PARPACK)

if (NOT @BUILD_SHARED_LIBS@)
	include(CMakeFindDependencyMacro)
	# Find dependencies
	if (NOT TARGET BLAS::BLAS)
		find_dependency(BLAS REQUIRED)
	endif()
	if (NOT TARGET LAPACK::LAPACK)
		find_dependency(LAPACK REQUIRED)
	endif()
	if (@ICB@)
		enable_language(Fortran)
	endif()
	if (@MPI@)
		include(FindMPI)
		if (NOT TARGET MPI::Fortran)
			find_dependency(MPI REQUIRED COMPONENTS Fortran)
		endif()
	endif()
endif()

include("${CMAKE_CURRENT_LIST_DIR}/arpackngTargets.cmake")

add_library(ARPACK::ARPACK ALIAS arpack)
if (TARGET parpack)
	add_library(PARPACK::PARPACK ALIAS parpack)
endif()
