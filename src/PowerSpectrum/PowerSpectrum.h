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
protected:
  GreensTensor *greens_tensor; // Green's tensor of describing the geometry of
  // the system
  Polarizability *polarizability; // Polarizability describing the linear
public:
  // response of the microscopic particle
  // Constructors
  PowerSpectrum(std::string input_file);

  // Constructor with initialization list
  PowerSpectrum(GreensTensor *greens_tensor, Polarizability *polarizability);

  // Calculate the power spectrum for a fixed value of the frequency \omega
  virtual void calculate(cx_mat::fixed<3, 3> &powerspectrum,
                         Options_PowerSpectrum opts) = 0;

  // getter functions
  GreensTensor *get_greens_tensor() { return greens_tensor; };
  Polarizability *get_polarizability() { return polarizability; };
};

struct Options_PowerSpectrum {
  // Compute the complete spectrum with LTE and non-LTE contributions
  bool full_spectrum = false;

  // Only compute the non-LTE contributions to the powerspectrum
  bool non_LTE = false;

  double omega = NAN;
};

#endif // POWERSPECTRUM_H
