#ifndef MEMORYKERNELFACTORY_H
#define MEMORYKERNELFACTORY_H

#include "MemoryKernel.h"

//! A Memory Kernel factory
/*!
* This is a class implementing the factory design pattern for memory kernels.
* Given an input file it returns a pointer to the right memory kernel.
* Possible options include: ohmic.
*/
class MemoryKernelFactory
{
public:
  /*!
  * Function returning a memory kernel pointer of the right type.
  * @param type Type of the memory kernel.
  */
  static MemoryKernel * create(std::string type);
};


#endif // MEMORYKERNELFACTORY_H
