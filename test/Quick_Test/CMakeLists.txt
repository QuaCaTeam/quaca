# add sources to test
file(GLOB_RECURSE test_sources ${CMAKE_CURRENT_SOURCE_DIR} "*.cpp")

# Executable
add_executable(quick_test
        quick_test.cpp
        ../UnitTests/Polarizability/test_PolarizabilityNoBathMagnetic_unit.cpp
	../UnitTests/GreensTensor/test_GreensTensorPlateMagnetic_unit.cpp
	../UnitTests/PowerSpectrum/test_PowerSpectrumHarmOscMagnetic.cpp
        )

target_include_directories(quick_test PRIVATE ../include)

target_link_libraries(quick_test PRIVATE
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
