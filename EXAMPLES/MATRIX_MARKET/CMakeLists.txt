if (ICB)
	if (ICBEXMM)
		add_arpack_test(
			USE_BASH
			NAME arpack_mm
			SOURCES arpackmm.cpp
			DATA_FILES
			An.mtx Az.mtx B.mtx Bz.mtx arpackmm.sh
			COMMAND arpackmm.sh)
	endif()
endif()
