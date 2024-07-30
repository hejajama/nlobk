configure_file("${CMAKE_CURRENT_SOURCE_DIR}/gitsha1.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/gitsha1.cpp" @ONLY)


add_executable(
	nlobk	
	main.cpp
	dipole.cpp
	solver.cpp
	mv.cpp
	ic.cpp
	ic_datafile.cpp
	nlobk_config.cpp
	../tools/tools.cpp
	../tools/interpolation.cpp
	${CMAKE_CURRENT_BINARY_DIR}/gitsha1.cpp
	)

include_directories(.)

find_package(GSL REQUIRED)    #
target_link_libraries(
	nlobk
	PRIVATE
#	amplitude
	GSL::gsl
	GSL::gslcblas
)


