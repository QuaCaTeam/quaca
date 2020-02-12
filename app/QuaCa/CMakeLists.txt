# add executable
add_executable(QuaCa
  ${CMAKE_CURRENT_SOURCE_DIR}/quaca.cpp
  )
target_include_directories(Friction PRIVATE ../include)

# link libraries
target_link_libraries(Friction
  quaca
  )
