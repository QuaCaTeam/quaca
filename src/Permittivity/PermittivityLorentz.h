#ifndef PERMITTIVITYLORENTZ_H
#define PERMITTIVITYLORENTZ_H

#include "../MemoryKernel/MemoryKernel.h"
#include "Permittivity.h"
#include <complex>

//! A Lorentz model permittivity
class PermittivityLorentz : public Permittivity {
private:
  double eps_inf;
  double omega_p;
  double omega_0;

  MemoryKernel *memory_kernel;

public:
  // constructors
  PermittivityLorentz(double eps_inf, double omega_p, double omega_0,
                      MemoryKernel *memory_kernel);
  PermittivityLorentz(const std::string& input_file);

  // calculate the permittivity
  std::complex<double> epsilon(double omega);

  // Returns the numerical value of the permittivity scaled by omega.
  std::complex<double> epsilon_omega(double omega);

  // getter methods
  double get_eps_inf() const { return this->eps_inf; };
  double get_omega_p() const { return this->omega_p; };
  double get_omega_0() const { return this->omega_0; };
  MemoryKernel *get_memory_kernel() const { return this->memory_kernel; };
};

#endif // PERMITTIVITYLORENTZ_H
