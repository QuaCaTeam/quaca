#ifndef PERMITTIVITY_H
#define PERMITTIVITY_H

#include <complex>

//! An abstract permittivity class
/*!
* This is an abstract class for permittivities
* All permittivities should return a complex number, given a real frequency as input.
* Being susceptibilities, all permittivities should obey the crossing relation.
*/
class Permittivity
{
public:

  /*!
  * Return the permittivity given at frequency omega
  * @param omega Frequency
  */
  virtual std::complex<double> epsilon(double omega) =0;
};

#endif // PERMITTIVITY_H
