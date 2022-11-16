#IF lhapdf-config DOESN'T WORK FOR YOU, CHECK THE TWO COMMENTS BELOW
add_library(
	Parts ars.cpp
	hcubature.cpp
	histo.cpp
	linear_interpolator.cpp
	mcaa.cpp
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
	"${lhapdf_incdir}" #REPLACE THIS WITH THE PATH TO THE DIRECTORY OF LHAPDF HEADER FILES (.h)!!!
)

target_compile_features(Parts PUBLIC cxx_std_20)

target_compile_options(Parts PUBLIC
	-Wall
	#-fconcepts #COMMENT THIS IF YOU USE CLANG
)

find_package(Boost REQUIRED)
find_package(GSL REQUIRED)

target_link_directories(Parts PUBLIC 
	"${lhapdf_libdir}" #REPLACE THIS WITH THE PATH TO THE DIRECTORY OF LHAPDF LIBRARY FILES (.so)!!!
)

target_link_libraries(Parts PUBLIC 
	GSL::gsl
	GSL::gslcblas
	-ltbb
	-lpthread
	-lm
	-lLHAPDF
)
