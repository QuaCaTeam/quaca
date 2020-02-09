# add executable
add_executable(Tutorial
  ${CMAKE_CURRENT_SOURCE_DIR}/tutorial.cpp
  )
target_include_directories(Tutorial PRIVATE ../include)

# link libraries
target_link_libraries(Tutorial
  quaca
  )
