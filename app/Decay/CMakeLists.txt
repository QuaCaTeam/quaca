# add executable
add_executable(Decay
  decay.cpp
  )
target_include_directories(Decay PRIVATE ../include)

# link libraries
target_link_libraries(Decay
  quaca
  OpenMP::OpenMP_CXX
  )
#OpenMP::OpenMP_CXX

