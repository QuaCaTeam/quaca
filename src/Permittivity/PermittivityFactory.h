#ifndef PERMITTIVITYFACTORY_H
#define PERMITTIVITYFACTORY_H

#include "Permittivity.h"
#include <memory>

//! A Permittivity factory
class PermittivityFactory {
public:
  // Returns a permittivity pointer of the right type.
  static std::shared_ptr<Permittivity> create(const std::string& type);
};

#endif // PERMITTIVITYFACTORY_H
