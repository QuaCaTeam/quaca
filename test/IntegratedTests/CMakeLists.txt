# add sources to test
set(test_sources
        test_main.cpp
        Friction/test_Friction_integrated.cpp
        Polarizability/test_PolarizabilityNoBath_integrated.cpp
        Polarizability/test_PolarizabilityBath_integrated.cpp
        )

# Executable
add_executable(test_quaca_integrated
  ${test_sources}
  )

target_include_directories(test_quaca_integrated PRIVATE ../include)

target_link_libraries(test_quaca_integrated PRIVATE
  catch2
  quaca
  ${GSL_LIBRARY}
  ${GSL_CBALS_LIBRARY}
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
)
