if (MPI)
	add_arpack_examples(
		LIBRARIES ${parpack_target} MPI::MPI_Fortran
		SOURCES
		pcndrv1.f
		pdndrv1.f
		pdndrv3.f
		pdsdrv1.f
		psndrv1.f
		psndrv3.f
		pssdrv1.f
		pzndrv1.f)
endif()
