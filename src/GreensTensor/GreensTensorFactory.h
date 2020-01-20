#ifndef GREENSTENSORFACTORY_H
#define GREENSTENSORFACTORY_H

#include "GreensTensor.h"

//! A Greens tensor factory
class GreensTensorFactory {
public:
  // Function returning a memory kernel pointer of the right type.
  static GreensTensor *create(std::string type);
};

#endif // GREENSTENSORFACTORY_H
