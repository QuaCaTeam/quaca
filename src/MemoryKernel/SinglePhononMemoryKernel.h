#ifndef SINGLEPHONONMEMORYKERNEL_H
#define SINGLEPHONONMEMORYKERNEL_H

#include "MemoryKernel.h"
#include <cmath>
#include <complex>

//! An SinglePhonon-shape memory kernel
class SinglePhononMemoryKernel : public MemoryKernel {
private:
  double gamma;        // ohmic damping coefficient
  double gamma_phon;   // phononic damping coefficient
  double omega_phon;   // phononic frequency
  double coupling;     // coupling coefficient to the dipole moment

public:
  // direct constructor
  explicit SinglePhononMemoryKernel(double gamma, double gamma_phon, double omega_phon, double coupling);

  // constructor from .json file
  explicit SinglePhononMemoryKernel(const std::string &input_file);

  // constructor from .json file of a specific section
  SinglePhononMemoryKernel(const std::string &input_file, const std::string &section);


  // calculate function
  std::complex<double> calculate(double omega) const override;

  // getter functions
  double get_gamma() const { return this->gamma; };
  double get_gamma_phon() const { return this->gamma_phon; };
  double get_omega_phon() const { return this->omega_phon; };
  double get_coupling() const { return this->coupling; };
};

#endif // SINGLEPHONONMEMORYKERNEL_H
