# add executable
add_executable(Friction
  ${CMAKE_CURRENT_SOURCE_DIR}/friction.cpp
  )
target_include_directories(Friction PRIVATE ../include)

# link libraries
target_link_libraries(Friction
  quaca
  )
