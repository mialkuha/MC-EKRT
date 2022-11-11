cmake_minimum_required(VERSION 3.10)

project(mcaa CXX)

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

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")

add_subdirectory(src)

add_executable(mcaa main.cpp)

target_compile_features(mcaa PUBLIC cxx_std_20)

target_compile_options(mcaa PUBLIC
	-Wall
	-fconcepts
)

target_include_directories(mcaa PUBLIC 
	"${PROJECT_BINARY_DIR}"
	"${lhapdf_incdir}"
)

find_package(Boost REQUIRED)
find_package(GSL REQUIRED)

target_link_directories(mcaa PUBLIC 
	"${lhapdf_libdir}"
	"/usr/local/lib" #TODO INCLUDE BOOST NORMALLY
)

target_link_libraries(mcaa PUBLIC 
	GSL::gsl
	GSL::gslcblas
	-ltbb
	-lpthread
	-lm
	-lLHAPDF
	Parts
)