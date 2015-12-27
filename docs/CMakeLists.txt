# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/TraceHBonds.in
		${CMAKE_CURRENT_BINARY_DIR}/TraceHBonds.cfg @ONLY)
	add_custom_target(docs
		${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/TraceHBonds.cfg
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
		COMMENT "Generating API documentation with Doxygen" VERBATIM
	)

	FILE(GLOB image_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/.. ${CMAKE_CURRENT_SOURCE_DIR}/../images/*.png)

	# Copy necessary doc files
	FOREACH (file ${image_files})
		CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/../${file} ${CMAKE_CURRENT_BINARY_DIR}/html/${file} COPYONLY)
	ENDFOREACH (file)
endif(DOXYGEN_FOUND)




