#ifndef POWERSPECTRUM_H
#define POWERSPECTRUM_H

#include "../GreensTensor/GreensTensor.h"
#include "../Polarizability/Polarizability.h"

struct Options_PowerSpectrum;
/*!
 *  This is an abstract class implementing a structure to compute the power
 * spectrum tensor
 */
class PowerSpectrum {
public:
    
  //Green's tensor of describing the geometry of
  //the system
  GreensTensor *greens_tensor; 

  //Polarizability describing the linear of the microscopic particle
  Polarizability *polarizability; 

  //Constructor with ini-file
  PowerSpectrum(std::string input_file);

  //Constructor with initialization list
  PowerSpectrum(GreensTensor *greens_tensor, Polarizability *polarizability);

  //Calculate the power spectrum for a fixed value of the frequency \omega
  virtual void calculate(cx_mat::fixed<3, 3> &powerspectrum,
                         Options_PowerSpectrum opts) = 0;
};

struct Options_PowerSpectrum {
  //Compute the complete spectrum with LTE and non-LTE contributions
  bool full_spectrum = false;

  //Only compute the non-LTE contributions to the powerspectrum
  bool non_LTE = false;

  double omega = NAN;
};

#endif // POWERSPECTRUM_H
