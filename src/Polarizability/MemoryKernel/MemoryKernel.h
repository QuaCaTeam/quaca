#ifndef MEMORYKERNEL_H
#define MEMORYKERNEL_H

#include <complex>
#include <cmath>

//! An abstract class for memory kernels
/*!
* This is an abstract class for memory kernels.
* All memory kernels should return a complex number, given a real frequency as input.
* Furthermore, all memory kernels should obey the crossing relation.
*/
class MemoryKernel
{
public:

    /*!
    * Returns the memory kernel given a frequency omega.
    * @param omega Frequency
    */
    virtual std::complex<double> mu( double omega ) =0;
};

#endif //MEMORYKERNEL_H
