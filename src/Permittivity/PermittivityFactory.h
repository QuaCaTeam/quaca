#ifndef PERMITTIVITYFACTORY_H
#define PERMITTIVITYFACTORY_H

#include <memory>
#include "Permittivity.h"

//! A Permittivity factory
class PermittivityFactory {
public:
  // Returns a permittivity pointer of the right type.
  static std::shared_ptr<Permittivity> create(const std::string &input_file);
};

#endif // PERMITTIVITYFACTORY_H
