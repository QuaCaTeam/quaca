#ifndef POWERSPECTRUM_H
#define POWERSPECTRUM_H

#include "../GreensTensor/GreensTensor.h"
#include "../Polarizability/Polarizability.h"

struct Options_PowerSpectrum;
enum Spectrum_Options { FULL, NON_LTE_ONLY };
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
  PowerSpectrum(const std::string &input_file);

  // Constructor with initialization list
  PowerSpectrum(GreensTensor *greens_tensor, Polarizability *polarizability);

  // Calculate the power spectrum for a fixed value of the frequency \omega
  virtual void calculate(cx_mat::fixed<3, 3> &powerspectrum,
                         Options_PowerSpectrum opts) = 0;

  // getter functions
  GreensTensor *get_greens_tensor() const { return greens_tensor; };
  Polarizability *get_polarizability() const { return polarizability; };
};

struct Options_PowerSpectrum {
  // Compute the complete spectrum with LTE and non-LTE contributions
  Spectrum_Options spectrum = FULL;

  // Only compute the non-LTE contributions to the powerspectrum

  double omega = NAN;
};

#endif // POWERSPECTRUM_H
