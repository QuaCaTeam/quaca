# add sources to test
file(GLOB_RECURSE test_sources ${CMAKE_CURRENT_SOURCE_DIR} "*.cpp")

# Executable
add_executable(test_quaca
  ${test_sources}
  )

target_include_directories(test_quaca PRIVATE ../include)

target_link_libraries(test_quaca PRIVATE
  catch2
  quaca
  gsl
  )
