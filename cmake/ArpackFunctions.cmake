include_guard()

find_program(BASH_PROGRAM bash)

function(add_arpack_test)
	#
	# add_arpack_test(
	#   [USE_BASH|USE_MPI]
	#   [NAME test_name]
	#   [OUTPUT_DIRECTORY path_in_build_directory]
	#   [COMMAND executable_or_target_name]
	#   [SOURCES source1 ... ]
	#   [LIBRARIES library1 ... ]
	#   [DATA_FILES file1 ... ])
	#
	# Fabricate a test with possible input data files, copying
	# them to the build directory, and using the right sources
	# and libraries.
	#
	cmake_parse_arguments(
		TEST
		"USE_BASH;USE_MPI" # options
		"NAME;OUTPUT_DIRECTORY;COMMAND" # 1 arg variables
		"SOURCES;LIBRARIES;DATA_FILES" # multiargument
		${ARGN})
	if ("${TEST_NAME}" STREQUAL "")
		list(GET TEST_SOURCES 0 first_source)
		get_filename_component(TEST_NAME ${first_source} NAME_WE)
		set(executable_name ${TEST_NAME}_test)
	else()
		set(executable_name ${TEST_NAME})
	endif()
	set(test_name ${TEST_NAME}_tst)

	if ("${TEST_COMMAND}" STREQUAL "")
		set(TEST_COMMAND ${executable_name})
	endif()
	if (MINGW)
		# Unfortunately, CMake in MINGW has problems with
		# the location of the library
		set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
	endif()

	add_executable(${executable_name} ${TEST_SOURCES})
	target_link_libraries(${executable_name} ${TEST_LIBRARIES} ${arpack_target})

	# Copy data files to the location where the test is built and run
	foreach(f ${TEST_DATA_FILES})
		file(COPY "${f}" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
	endforeach()

	if (TEST_USE_BASH)
		add_test(NAME ${test_name}
			COMMAND ${BASH_PROGRAM} ${TEST_COMMAND}
			WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
	else()
		if (TEST_USE_MPI)
			add_test(NAME ${test_name}
				COMMAND mpirun -n 2 ${TEST_COMMAND}
				WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
		else()
			add_test(NAME ${test_name}
				COMMAND ${executable_name}
				WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
		endif()
	endif()
endfunction()

function(add_pyarpack_test template_file_name)
	#
	# make_pyarpack_test(full_path_to_py_in_file)
	#
	# convert a *.py.in file to Pyhon and create a test that
	# executes this script, including pyarpack in the path.
	#
	get_filename_component(name_root ${template_file_name} NAME_WE)
	set(test_name "${name_root}_tst")
	set(script_name "${CMAKE_BINARY_DIR}/TESTS/${name_root}.py")
	configure_file("${template_file_name}" "${script_name}" @ONLY)
	add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} ${script_name})
	set_tests_properties(${test_name}
		PROPERTIES ENVIRONMENT
		PYTHONPATH=${CMAKE_BINARY_DIR}/lib:$ENV{PYTHONPATH}
	)
endfunction()

function(add_arpack_examples)
	#
	# add_arpack_examples(
	#    [SOURCES test1 test2 ...]
	#    [COMMON_SOURCES file1 file2 ...]
	#    [LIBRARIES lib1 lib2 ...]
	# )
	#
	# Create a set of tests for each of the files in SOURCES,
	# possibly adding LIBRARIES and other shared COMMON_SOURCES
	#
	cmake_parse_arguments(
		EXAMPLE
		"" # options
		"COMMON_SOURCES" # 1 arg variables
		"SOURCES;LIBRARIES" # multiargument
		${ARGN})
	foreach(file ${EXAMPLE_SOURCES})
		get_filename_component(example_name ${file} NAME_WE)
		add_arpack_test(NAME ${example_name}
			SOURCES "${file};${EXAMPLE_COMMON_SOURCES}"
			LIBRARIES "${EXAMPLE_LIBRARIES}"
		)
	endforeach()
endfunction()
