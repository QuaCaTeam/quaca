#ifndef GREENSTENSORFACTORY_H
#define GREENSTENSORFACTORY_H

#include "GreensTensor.h"
#include <memory>

//! A Greens tensor factory
/*!
 * This is a class implementing the factory design pattern for Green's tensors.
 * Given an input file it returns a pointer to the right Green's tensor.
 * Possible options include: vacuum, plate.
 */
class GreensTensorFactory {
public:
  /*!
   * Function returning a memory kernel pointer of the right type.
   * @param type Type of the memory kernel.
   */
  static std::shared_ptr<GreensTensor> create(const std::string &input_file);
};

#endif // GREENSTENSORFACTORY_H
