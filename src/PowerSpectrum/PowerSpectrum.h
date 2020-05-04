#ifndef POWERSPECTRUM_H
#define POWERSPECTRUM_H

#include "../GreensTensor/GreensTensor.h"
#include "../Polarizability/Polarizability.h"

enum Spectrum_Options { FULL, NON_LTE_ONLY };

/*!
 *  This is an abstract class implementing a structure to compute the power
 * spectrum tensor
 */
class PowerSpectrum {
protected:
  std::shared_ptr<GreensTensor>
      greens_tensor; // Green's tensor of describing the geometry of
  // the system
  std::shared_ptr<Polarizability>
      polarizability; // Polarizability describing the linear
public:
  // response of the microscopic particle
  // Constructors
  PowerSpectrum(const std::string &input_file);

  // Constructor with initialization list
  PowerSpectrum(std::shared_ptr<GreensTensor> greens_tensor,
                std::shared_ptr<Polarizability> polarizability);

  // Calculate the power spectrum for a fixed value of the frequency \omega
  void calculate(double omega, cx_mat::fixed<3, 3> &powerspectrum,
                 Spectrum_Options spectrum) const;

  // getter functions
  std::shared_ptr<GreensTensor> &get_greens_tensor() { return greens_tensor; };
  std::shared_ptr<Polarizability> &get_polarizability() {
    return polarizability;
  };
};

#endif // POWERSPECTRUM_H
