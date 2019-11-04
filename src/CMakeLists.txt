# add QuaCa library
file(GLOB_RECURSE quaca_sources ${CMAKE_CURRENT_SOURCE_DIR} "*.cpp")
add_library(quaca SHARED ${quaca_sources})

# link boost libraries
target_link_libraries(quaca
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
)

# set global header
set_target_properties(quaca PROPERTIES PUBLIC_HEADER include/Quaca.h)
