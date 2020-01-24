# add QuaCa library
file(GLOB_RECURSE quaca_sources ${CMAKE_CURRENT_SOURCE_DIR} "*.cpp")
add_library(quaca SHARED ${quaca_sources})

# link libraries
target_link_libraries(quaca
    ${GSL_LIBRARY}
    ${GSL_CBLAS_LIBRARY}
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
)

# set global header
set_target_properties(quaca PROPERTIES PUBLIC_HEADER include/Quaca.h)
