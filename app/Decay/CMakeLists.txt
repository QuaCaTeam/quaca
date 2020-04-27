# add executable
add_executable(Decay
  ${CMAKE_CURRENT_SOURCE_DIR}/decay.cpp
  )
target_include_directories(Decay PRIVATE ../include)

# link libraries
target_link_libraries(Decay
  quaca
  OpenMP::OpenMP_CXX
  )
#OpenMP::OpenMP_CXX

