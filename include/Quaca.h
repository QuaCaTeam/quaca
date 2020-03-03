#ifndef QUACA_H
#define QUACA_H

#include "../src/Calculations/Integrations.h"

#include "../src/GreensTensor/GreensTensor.h"
#include "../src/GreensTensor/GreensTensorFactory.h"
#include "../src/GreensTensor/GreensTensorPlate.h"
#include "../src/GreensTensor/GreensTensorPlateMagnetic.h"
#include "../src/GreensTensor/GreensTensorPlateVacuum.h"
#include "../src/GreensTensor/GreensTensorVacuum.h"

#include "../src/Looper/Looper.h"
#include "../src/Looper/LooperFactory.h"
#include "../src/Looper/LooperV.h"
#include "../src/Looper/LooperZa.h"

#include "../src/MemoryKernel/MemoryKernel.h"
#include "../src/MemoryKernel/OhmicMemoryKernel.h"

#include "../src/Options/Options.h"

#include "../src/Permittivity/Permittivity.h"
#include "../src/Permittivity/PermittivityDrude.h"
#include "../src/Permittivity/PermittivityFactory.h"

#include "../src/Polarizability/Polarizability.h"
#include "../src/Polarizability/PolarizabilityFactory.h"
#include "../src/Polarizability/PolarizabilityBath.h"
#include "../src/Polarizability/PolarizabilityNoBath.h"
#include "../src/Polarizability/PolarizabilityNoBathMagnetic.h"

#include "../src/PowerSpectrum/PowerSpectrum.h"
#include "../src/PowerSpectrum/PowerSpectrumHarmOsc.h"
#include "../src/PowerSpectrum/PowerSpectrumHarmOscMagnetic.h"
#include "../src/PowerSpectrum/PowerSpectrumFactory.h"

#include "../src/Friction/Friction.h"

#include "../src/ReflectionCoefficients/ReflectionCoefficients.h"
#include "../src/ReflectionCoefficients/ReflectionCoefficientsFactory.h"
#include "../src/ReflectionCoefficients/ReflectionCoefficientsLocBulk.h"
#include "../src/ReflectionCoefficients/ReflectionCoefficientsLocSlab.h"

#endif //QUACA_H
