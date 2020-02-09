#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <iostream>
namespace pt = boost::property_tree;

#include "ProgressBar.hpp"
#include "Quaca.h"

int main(int argc, char *argv[]) {

  PermittivityDrude permittivity(9.0, 0.1);
  ReflectionCoefficientsLocBulk refl_coefficients(&permittivity);
  GreensTensorPlate greens_tensor(1e-4, 0.01, 1e6, &refl_coefficients, 20,
                                  {1e-4, 1e-2});
  PolarizabilityNoBath polarizability(1.3, 6e-9, &greens_tensor);
  PowerSpectrumHarmOsc power_spectrum(&greens_tensor, &polarizability);
  QuantumFriction friction(&greens_tensor, &polarizability, &power_spectrum,
                           1e-1);

  // quantum friction options
  Options_Friction opts;
  opts.non_LTE = true;
  opts.class_pt = &friction;

  // loop over v
  double start = 1e-4;
  double end = 1e-2;
  int number_of_steps = 2;
  double spacing = pow(end / start, 1. / ((double)number_of_steps - 1.0));
  double step, value;

  // define output file
  std::ofstream file;
  file.open("tutorial_mainfile.csv");

  for (int i = 0; i < number_of_steps; i++) {
    step = start * pow(spacing, i);
    friction.greens_tensor->set_v(step);
    value = friction.calculate(opts);

    std::cout << step << "," << value << std::endl;
    file << step << "," << value << "\n";
  };

  // close file
  file.close();

  return 0;
};
