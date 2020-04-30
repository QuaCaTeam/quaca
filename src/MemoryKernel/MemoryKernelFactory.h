#ifndef MEMORYKERNELFACTORY_H
#define MEMORYKERNELFACTORY_H

#include "MemoryKernel.h"
#include <memory>

//! A Memory Kernel factory
class MemoryKernelFactory {
public:
  // Function returning a memory kernel pointer of the right type.
  static std::shared_ptr<MemoryKernel> create(const std::string &input_file,
                                              const std::string &section);
};

#endif // MEMORYKERNELFACTORY_H
