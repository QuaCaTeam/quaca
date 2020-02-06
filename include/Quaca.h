#ifndef QUACA_H
#define QUACA_H

#include "../src/Calculations/Integrations.h"
#include "../src/Calculations/QuantumFriction.h"

#include "../src/GreensTensor/GreensTensor.h"
#include "../src/GreensTensor/GreensTensorFactory.h"
#include "../src/GreensTensor/GreensTensorPlate.h"
#include "../src/GreensTensor/GreensTensorVacuum.h"
#include "../src/GreensTensor/Permittivity/Permittivity.h"
#include "../src/GreensTensor/Permittivity/PermittivityDrude.h"
#include "../src/GreensTensor/Permittivity/PermittivityFactory.h"
#include "../src/GreensTensor/ReflectionCoefficients/ReflectionCoefficients.h"
#include "../src/GreensTensor/ReflectionCoefficients/ReflectionCoefficientsFactory.h"
#include "../src/GreensTensor/ReflectionCoefficients/ReflectionCoefficientsLocBulk.h"
#include "../src/GreensTensor/ReflectionCoefficients/ReflectionCoefficientsLocSlab.h"
#include "../src/Looper/Looper.h"
#include "../src/Looper/LooperFactory.h"
#include "../src/Looper/LooperV.h"
#include "../src/Looper/LooperZa.h"
#include "../src/Options/Options.h"
#include "../src/Polarizability/Polarizability.h"
#include "../src/Polarizability/PolarizabilityFactory.h"
#include "../src/Polarizability/PolarizabilityBath.h"
#include "../src/Polarizability/PolarizabilityNoBath.h"
#include "../src/Polarizability/MemoryKernel/MemoryKernel.h"
#include "../src/Polarizability/MemoryKernel/OhmicMemoryKernel.h"

#include "../src/PowerSpectrum/PowerSpectrum.h"
#include "../src/PowerSpectrum/PowerSpectrumHarmOsc.h"
#include "../src/PowerSpectrum/PowerSpectrumFactory.h"

#endif //QUACA_H
