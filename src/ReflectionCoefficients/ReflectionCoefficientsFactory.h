#ifndef REFLECTIONCOEFFICIENTSFACTORY_H
#define REFLECTIONCOEFFICIENTSFACTORY_H

#include <memory>
#include "ReflectionCoefficients.h"

//! A factory for the reflection coefficients
class ReflectionCoefficientsFactory {
public:
  /*!
   * Function returning a reflection coefficients pointer of the right type.
   * @param type Type of the reflection coefficients.
   */
  static std::shared_ptr<ReflectionCoefficients>
  create(const std::string &input_file);
};

#endif // REFLECTIONCOEFFICIENTSFACTORY_H
