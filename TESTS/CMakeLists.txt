add_arpack_test(SOURCES dnsimp.f mmio.f debug.h
	DATA_FILES testA.mtx
)

if (ICB)
    add_arpack_test(SOURCES bug_1315_single.c)
	add_arpack_test(SOURCES bug_1315_double.c)
endif()

add_arpack_test(SOURCES bug_1323.f)
add_arpack_test(SOURCES bug_58_double.f)
add_arpack_test(SOURCES bug_79_double_complex.f)
add_arpack_test(SOURCES bug_142.f)
add_arpack_test(SOURCES bug_142_gen.f)

if(ICB)
	add_arpack_test(SOURCES icb_arpack_c.c)
	add_arpack_test(SOURCES icb_arpack_cpp.cpp)
endif()
