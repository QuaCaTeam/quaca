#ifndef OHMICMEMORYKERNEL_H
#define OHMICMEMORYKERNEL_H

#include <complex>
#include <cmath>
#include "MemoryKernel.h"

//! An ohmic memory kernel
/*!
* This is a class implementing an ohmic memory kernel, i.e. a memory kernel
* given by \f$ \mu(\omega) = \gamma \f$.
* It takes the damping coefficient as an input
*/
class OhmicMemoryKernel : public MemoryKernel
{
  private:
    /*! damping coefficient */
    double gamma;

  public:

    /*!
    * Returns the memory kernel given a frequency omega, given by the formula
    * \f$ \mu(\omega) = \gamma \f$
    * @param omega Frequency
    */
    std::complex<double> mu( double omega );

    /*!
    * Constructor for an ohmic memory kernel
    * @param a damping coefficient
    */
    OhmicMemoryKernel(double a);

    /*!
    * Constructor for an ohmic memory kernel
    * @param input_file .ini file containing parameters
    */
    OhmicMemoryKernel(std::string input_file);

    /*!
    * Getter method for damping coefficient.
    */
    double get_gamma();

};

#endif //OHMICMEMORYKERNEL_H
