#ifndef MEMORYKERNELFACTORY_H
#define MEMORYKERNELFACTORY_H

#include "MemoryKernel.h"

//! A Memory Kernel factory
class MemoryKernelFactory {
public:
  // Function returning a memory kernel pointer of the right type.
  static MemoryKernel *create(std::string type);
};

#endif // MEMORYKERNELFACTORY_H
