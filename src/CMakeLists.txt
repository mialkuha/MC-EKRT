add_library(
	Parts ars.cpp
	hcubature.cpp
	histo.cpp
	linear_interpolator.cpp
	nn_coll.cpp
	nucleus_generator.cpp
	pqcd.cpp
	eps09.cxx
)

execute_process(
	COMMAND bash -c "lhapdf-config --incdir"
	OUTPUT_VARIABLE lhapdf_incdir
)
string(STRIP ${lhapdf_incdir} lhapdf_incdir)

execute_process(
	COMMAND bash -c "lhapdf-config --libdir"
	OUTPUT_VARIABLE lhapdf_libdir
) 
string(STRIP ${lhapdf_libdir} lhapdf_libdir)

target_include_directories(Parts PUBLIC 
	"${CMAKE_CURRENT_SOURCE_DIR}"
	"${lhapdf_incdir}"
)

target_compile_features(Parts PUBLIC cxx_std_20)

target_compile_options(Parts PUBLIC
	-Wall
	-fconcepts
)

find_package(Boost REQUIRED)
find_package(GSL REQUIRED)

target_link_directories(Parts PUBLIC 
	"${lhapdf_libdir}"
	"/usr/local/lib" #TODO INCLUDE BOOST NORMALLY
)

target_link_libraries(Parts PUBLIC 
	GSL::gsl
	GSL::gslcblas
	-ltbb
	-lpthread
	-lm
	-lLHAPDF
)
