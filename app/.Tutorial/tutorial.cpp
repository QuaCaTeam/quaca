#include "ProgressBar.hpp"
#include "Quaca.h"

int main(int argc, char *argv[]) {

  // parameters for permittivity
  double omega_p = 9.0;
  double gamma = 0.1;

  // define permittivity
  PermittivityDrude permittivity(omega_p, gamma);

  // define reflection coefficients
  ReflectionCoefficientsLocBulk refl_coefficients(&permittivity);

  // parameters for green's tensor
  double v = 1e-4;
  double beta = 1e6;
  double z_a = 0.01;

  // numerical error for green's tensor
  double delta_cut = 20;
  vec::fixed<2> rel_err = {1E-4, 1E-2};

  // define the Green's tensor
  GreensTensorPlate greens_tensor(v, z_a, beta, &refl_coefficients, delta_cut,
                                  rel_err);

  // parameters for polarizability
  double omega_a = 1.3;
  double alpha_zero = 6e-9;

  // define polarizability
  PolarizabilityNoBath polarizability(omega_a, alpha_zero, &greens_tensor);

  // define power spectrum
  PowerSpectrumHarmOsc power_spectrum(&greens_tensor, &polarizability);

  // numerical error for quantum friction
  double rel_err_omega = 1e-1;

  // define quantum friction
  Friction friction(&greens_tensor, &polarizability, &power_spectrum,
                    rel_err_omega);

  // quantum friction options
  Options_Friction opts;
  opts.spectrum = NON_LTE_ONLY;
  opts.class_pt = &friction;

  // loop over v
  double start = 1e-4;
  double end = 1e-2;
  int number_of_steps = 40;
  double spacing = pow(end / start, 1. / ((double)number_of_steps - 1.0));

  // define output file
  std::ofstream file;
  file.open("tutorial_mainfile.csv");

  double step, value;
  for (int i = 0; i < number_of_steps; i++) {
    step = start * pow(spacing, i);
    friction.get_greens_tensor()->set_v(step);
    value = friction.calculate(opts);

    file << step << "," << value << "\n";
  };

  // close file
  file.close();

  return 0;
};
