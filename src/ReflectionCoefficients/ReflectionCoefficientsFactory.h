#ifndef REFLECTIONCOEFFICIENTSFACTORY_H
#define REFLECTIONCOEFFICIENTSFACTORY_H

#include "ReflectionCoefficients.h"

//! A factory for the reflection coefficients
class ReflectionCoefficientsFactory {
public:
  /*!
   * Function returning a reflection coefficients pointer of the right type.
   * @param type Type of the reflection coefficients.
   */
  static ReflectionCoefficients *create(const std::string &type);
};

#endif // REFLECTIONCOEFFICIENTSFACTORY_H
