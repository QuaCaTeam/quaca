#ifndef MEMORYKERNELFACTORY_H
#define MEMORYKERNELFACTORY_H

#include "MemoryKernel.h"

//! A Memory Kernel factory
class MemoryKernelFactory {
public:
  // Function returning a memory kernel pointer of the right type.
  static MemoryKernel *create(const std::string &type,
                              const std::string &section);
};

#endif // MEMORYKERNELFACTORY_H
