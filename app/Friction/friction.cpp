#include "Quaca.h"
#include <iostream>

int main(int argc, char *argv[]) {
  // get file from command line
  Options opts(argc, argv);
  std::string parameters = opts.get_parameter_file();

  // define needed quantities
  GreensTensor *greens_tensor = GreensTensorFactory::create(parameters);
  Polarizability *polarizabilty = PolarizabilityFactory::create(parameters);
  PowerSpectrumHarmOsc *powerspectrum =
      new PowerSpectrumHarmOsc(greens_tensor, polarizabilty);
  QuantumFriction *quant_friction =
      new QuantumFriction(greens_tensor, polarizabilty, powerspectrum);

  // define looper
  Looper *looper = LooperFactory::create(parameters, quant_friction);

  // define output file
  std::ofstream file;
  file.open(opts.get_output_file());

  // calculate values
  for (int i = 0; i < looper->get_steps_total(); i++) {
    file << looper->get_step(i) << "," << looper->calculate_value(i) << "\n";
  };

  // close file
  file.close();

  // delete dynamically allocated memory
  delete[] greens_tensor;
  delete[] polarizabilty;
  delete[] powerspectrum;
  delete[] quant_friction;

  return 0;
};
