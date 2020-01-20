#ifndef POLARIZABILITYFACTORY_H
#define POLARIZABILITYFACTORY_H

#include "Polarizability.h"

//! A Greens tensor factory
class PolarizabilityFactory {
public:
  // Function returning a memory kernel pointer of the right type.
  static Polarizability *create(std::string type);
};

#endif // GREENSTENSORFACTORY_H
