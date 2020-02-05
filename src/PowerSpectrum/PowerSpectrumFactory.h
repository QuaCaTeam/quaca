#ifndef POWERSPECTRUMFACTORY_H
#define POWERSPECTRUMFACTORY_H

#include "PowerSpectrum.h"

//! A Greens tensor factory
/*!
 * This is a class implementing the factory design pattern for the power
 * spectrum. Given an input file it returns a pointer to the right power
 * spectrum. Possible options include: harmonic oscillator.
 */
class PowerSpectrumFactory {
public:
  /*!
   * Function returning a memory kernel pointer of the right type.
   * @param type Type of the memory kernel.
   */
  static PowerSpectrum *create(std::string type);
};

#endif // POWERSPECTRUMFACTORY_H
