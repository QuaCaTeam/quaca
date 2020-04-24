cmake_minimum_required(VERSION 3.7.2)
project(QuaCa VERSION 1.0.0 LANGUAGES CXX)

# compile flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1z -Ofast")

# define output directories for library files and binaries
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# find boost libraries
set(Boost_NO_BOOST_CMAKE ON)
find_package(Boost COMPONENTS program_options filesystem system REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# find gsl
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

# find lapack and use armadillo without runtimelibrary, wrap blas and lapack
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)
include_directories("${PROJECT_SOURCE_DIR}/lib/armadillo-9.870.2/include")

#find OpenMP library
find_package(OpenMP REQUIRED)
#add -fopenmp to the cxx_flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DARMA_DONT_USE_WRAPPER")

# set include directories
include_directories(
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_BINARY_DIR}/include
  )

# catch2
add_library(catch2 INTERFACE)
target_include_directories(catch2 INTERFACE ${PROJECT_SOURCE_DIR}/include)

# Subdirectories
add_subdirectory(src)
add_subdirectory(app)
add_subdirectory(test)
