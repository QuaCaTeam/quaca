#ifndef POWERSPECTRUM
#define POWERSPECTRUM

#include "../GreensTensor/GreensTensor.h"
#include "../Polarizability/Polarizability.h"

/*!
 *  This is an abstract class implementing a structure tom compute the power spectrum tensor
 */

class PowerSpectrum
{
  protected:
    GreensTensor* greens_tensor; // Green's tensor of describing the geometry of the system
    Polarizability* polarizability; // Polarizability describing the linear response of the microscopic particle
  public:
    //Constructors
    PowerSpectrum(std::string input_file);
    PowerSpectrum(GreensTensor* greens_tensor, Polarizability* polarizability);

    // calculate the power spectrum for a fixed value of the frequency
    virtual void calculate(cx_mat::fixed<3,3>& powerspectrum, double omega)=0;
};

#endif
