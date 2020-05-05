#ifndef LORENTZMEMORYKERNEL_H
#define LORENTZMEMORYKERNEL_H

#include "MemoryKernel.h"
#include <cmath>
#include <complex>

//! An Lorentz-shape memory kernel
class LorentzMemoryKernel : public MemoryKernel {
private:
  double gamma;   // damping coefficient
  double omega_0; // central frequency
  double omega_p; // plasma frequency
  double eps_inf;  // high-frequency permittivity limit

public:
  // direct constructor
  explicit LorentzMemoryKernel(double gamma, double omega_0, double omega_p, double eps_inf);

  // constructor from .json file
  explicit LorentzMemoryKernel(const std::string &input_file);

  // constructor from .json file of a specific section
  LorentzMemoryKernel(const std::string &input_file, const std::string &section);


  // calculate function
  std::complex<double> calculate(double omega) const override;

  // getter functions
  double get_gamma() const { return this->gamma; };
  double get_omega_0() const { return this->omega_0; };
  double get_omega_p() const { return this->omega_p; };
  double get_eps_inf() const { return this->eps_inf; };
};

#endif // LORENTZMEMORYKERNEL_H
