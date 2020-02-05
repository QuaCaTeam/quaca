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
  std::cout << "break" << std::endl;

  // define output file
  std::ofstream file;
  file.open(opts.get_output_file());

  // calculate values
  double step, value;
  for (int i = 0; i < looper->get_steps_total(); i++) {
    step = looper->get_step(i);
    value = looper->calculate_value(i);

    std::cout << step << "," << value << std::endl;
    file << step << "," << value << "\n";
  };

  // close file
  file.close();

  // delete dynamically allocated memory
  // delete[] greens_tensor;
  // delete[] polarizabilty;
  // delete[] powerspectrum;
  // delete[] quant_friction;

  return 0;
};
