#ifndef POWERSPECTRUMFACTORY_H
#define POWERSPECTRUMFACTORY_H

#include "PowerSpectrum.h"

//! A Power spectrum factory
class PowerSpectrumFactory {
public:
  // Function returning a power spectrum pointer of the right type.
  static PowerSpectrum *create(std::string type);
};

#endif // POWERSPECTRUMFACTORY_H
