#ifndef PERMITTIVITYFACTORY_H
#define PERMITTIVITYFACTORY_H

#include "Permittivity.h"

//! A Permittivity factory
class PermittivityFactory {
public:
  // Returns a permittivity pointer of the right type.
  static Permittivity *create(std::string type);
};

#endif // PERMITTIVITYFACTORY_H
