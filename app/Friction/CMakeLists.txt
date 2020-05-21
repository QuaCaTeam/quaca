# add executable
add_executable(Friction
  friction.cpp
  )
target_include_directories(Friction PRIVATE ../include)

# link libraries
target_link_libraries(Friction
  quaca
  OpenMP::OpenMP_CXX
  )
#OpenMP::OpenMP_CXX

