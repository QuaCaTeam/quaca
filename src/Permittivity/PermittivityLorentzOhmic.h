#ifndef PERMITTIVITYLORENTZOHMIC_H
#define PERMITTIVITYLORENTZOHMIC_H

#include "Permittivity.h"
#include <complex>

//! A LorentzOhmic model permittivity
class PermittivityLorentzOhmic : public Permittivity {
private:
  double eps_inf;
  double alpha_0;
  double omega_0;
  double gamma;

public:
  // constructors
  PermittivityLorentzOhmic(double eps_inf, double alpha_zero, double omega_0,
                           double gamma);
  PermittivityLorentzOhmic(const std::string &input_file);

  // calculate the permittivity
  std::complex<double> calculate(double omega) const override;
  // calculate the permittivy times omega
  std::complex<double> calculate_times_omega(double omega) const override;

  // getter methods
  double get_eps_inf() const { return this->eps_inf; };
  double get_alpha_0() const { return this->alpha_0; };
  double get_omega_0() const { return this->omega_0; };
  double get_gamma() const { return this->gamma; };

  // print all physical quantities
  void print_info(std::ostream &stream) const override;
};
#endif
