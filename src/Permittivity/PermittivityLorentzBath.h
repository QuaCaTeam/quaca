#ifndef PERMITTIVITYLORENTZBATH_H
#define PERMITTIVITYLORENTZBATH_H

#include "../MemoryKernel/MemoryKernel.h"
#include "Permittivity.h"
#include <complex>

//! A LorentzBath model permittivity
class PermittivityLorentzBath : public Permittivity {
private:
  double eps_inf;
  double alpha_zero;
  double omega_0;

  MemoryKernel *memory_kernel;

public:
  // constructors
  PermittivityLorentzBath(double eps_inf, double alpha_zero, double omega_0,
                          MemoryKernel *memory_kernel);
  PermittivityLorentzBath(std::string input_file);

  // calculate the permittivity
  std::complex<double> epsilon(double omega);

  // Returns the numerical value of the permittivity scaled by omega.
  std::complex<double> epsilon_omega(double omega);

  // getter methods
  double get_eps_inf() { return this->eps_inf; };
  double get_alpha_zero() { return this->alpha_zero; };
  double get_omega_0() { return this->omega_0; };
  MemoryKernel *get_memory_kernel() { return this->memory_kernel; };
};

#endif // PERMITTIVITYLORENTZBATH_H
