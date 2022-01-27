if (MPI)
	add_arpack_test(SOURCES issue46.f
		LIBRARIES ${parpack_target} MPI::MPI_Fortran
	)

	if (ICB)
		add_arpack_test(
			USE_MPI
			SOURCES icb_parpack_c.c
			LIBRARIES ${parpack_target} MPI::MPI_C
		)
		if (ICBEXMM)
			add_arpack_test(
				USE_MPI
				SOURCES icb_parpack_cpp.cpp
				LIBRARIES ${parpack_target} MPI::MPI_CXX
			)
		endif()
	endif()
endif()
