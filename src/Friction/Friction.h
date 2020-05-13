#ifndef QUANTUMFRICTION_H
#define QUANTUMFRICTION_H

#include "../GreensTensor/GreensTensor.h"
#include "../Polarizability/Polarizability.h"
#include "../PowerSpectrum/PowerSpectrum.h"
#include <string>

/*!
 * This is a class computing the quantum friction force for a given Green's
 * tensor and polarizibility with a certain set of parameters
 */

class Friction {
protected:
  std::shared_ptr<GreensTensor> greens_tensor;
  std::shared_ptr<Polarizability> polarizability;
  std::shared_ptr<PowerSpectrum> powerspectrum;

  double relerr_omega;

public:
  Friction(const std::string &input_file);
  Friction(std::shared_ptr<GreensTensor> greens_tensor,
           std::shared_ptr<Polarizability> polarizability,
           std::shared_ptr<PowerSpectrum> powerspectrum, double relerr_omega);

  double calculate(Spectrum_Options spectrum) const;
  double friction_integrand(double omega, Spectrum_Options spectrum) const;

  // getter functions
  std::shared_ptr<GreensTensor> &get_greens_tensor() { return greens_tensor; };
  std::shared_ptr<Polarizability> &get_polarizability() {
    return polarizability;
  };
  std::shared_ptr<PowerSpectrum> &get_powerspectrum() { return powerspectrum; };
};

#endif // QUANTUMFRICTION_H
