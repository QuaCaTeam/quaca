#ifndef PERMITTIVITYFACTORY_H
#define PERMITTIVITYFACTORY_H

#include "Permittivity.h"

//! A Permittivity factory
/*!
* This is a class implementing the factory design pattern for permittivities.
* Given an input file it returns a pointer to the right permittivity.
* Possible options include: drude.
*/
class PermittivityFactory
{
public:
  /*!
  * Function returning a permittivity pointer of the right type.
  * @param type Type of the permittivity.
  */
  static Permittivity * create(std::string type);
};


#endif // PERMITTIVITYFACTORY_H
