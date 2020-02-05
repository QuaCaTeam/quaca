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
  GreensTensor *greens_tensor; // Green's tensor of describing the geometry of
  // the system
  Polarizability *polarizability; // Polarizability describing the linear
  // response of the microscopic particle
  // Constructors
  PowerSpectrum(std::string input_file);
  PowerSpectrum(GreensTensor *greens_tensor, Polarizability *polarizability);

  // calculate the power spectrum for a fixed value of the frequency
  virtual void calculate(cx_mat::fixed<3, 3> &powerspectrum,
                         Options_PowerSpectrum opts) = 0;
};

struct Options_PowerSpectrum {
  bool full_spectrum = false;
  bool non_LTE = false;
  double omega = NAN;
};

#endif // POWERSPECTRUM_H
