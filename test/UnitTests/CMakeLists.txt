# add sources to test
set(test_sources
        test_main.cpp
        Calculations/test_Integrations_unit.cpp
        Friction/test_Friction_unit.cpp
        GreensTensor/test_GreensTensorPlate_unit.cpp
        GreensTensor/test_GreensTensorPlateVacuum_unit.cpp
        GreensTensor/test_GreensTensorVacuum_unit.cpp
        Looper/test_LooperV_unit.cpp
        Looper/test_LooperZa_unit.cpp
        Looper/test_LooperBeta_unit.cpp
        MemoryKernel/test_OhmicMemoryKernel_unit.cpp
        MemoryKernel/test_SinglePhononMemoryKernel_unit.cpp
        Permittivity/test_PermittivityDrude_unit.cpp
        Permittivity/test_PermittivityLorentz_unit.cpp
        Polarizability/test_PolarizabilityBath_unit.cpp
        Polarizability/test_PolarizabilityNoBath_unit.cpp
        PowerSpectrum/test_PowerSpectrum.cpp
        ReflectionCoefficients/test_ReflectionCoefficientsLocBulk_unit.cpp
        ReflectionCoefficients/test_ReflectionCoefficientsLocSlab_unit.cpp
        )

# Executable
add_executable(test_quaca_unit
        ${test_sources}
        )

target_include_directories(test_quaca_unit PRIVATE ../include)

target_link_libraries(test_quaca_unit PRIVATE
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
