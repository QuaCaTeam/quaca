cmake_minimum_required(VERSION 3.7.2)
project(QuaCa VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1z -Ofast")

#set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR})
#set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# find boost package
set(Boost_NO_BOOST_CMAKE ON)
find_package(Boost COMPONENTS filesystem system REQUIRED)

# use armadillo without runtimelibrary, wrap blas and lapack
include_directories("${PROJECT_SOURCE_DIR}/lib/armadillo-9.800.2/include")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DARMA_DONT_USE_WRAPPER -lblas -llapack")

# include directories
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
