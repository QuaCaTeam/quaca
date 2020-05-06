# add QuaCa library
set(quaca_sources
        Calculations/Integrations.cpp
        Friction/Friction.cpp
        GreensTensor/GreensTensor.cpp
        GreensTensor/GreensTensorFactory.cpp
        GreensTensor/GreensTensorPlate.cpp
        GreensTensor/GreensTensorPlateVacuum.cpp
        GreensTensor/GreensTensorVacuum.cpp
        Looper/Looper.cpp
        Looper/LooperFactory.cpp
        Looper/LooperV.cpp
        Looper/LooperZa.cpp
        MemoryKernel/MemoryKernelFactory.cpp
        MemoryKernel/OhmicMemoryKernel.cpp
	MemoryKernel/SinglePhononMemoryKernel.cpp
        Options/Options.cpp
        Permittivity/PermittivityDrude.cpp
        Permittivity/PermittivityFactory.cpp
        Permittivity/PermittivityLorentz.cpp
        Polarizability/Polarizability.cpp
        PowerSpectrum/PowerSpectrum.cpp
        ReflectionCoefficients/ReflectionCoefficientsFactory.cpp
        ReflectionCoefficients/ReflectionCoefficientsLocBulk.cpp
        ReflectionCoefficients/ReflectionCoefficientsLocSlab.cpp
        )
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
