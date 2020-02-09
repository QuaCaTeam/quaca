#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <iostream>
namespace pt = boost::property_tree;

#include "ProgressBar.hpp"
#include "Quaca.h"

int main(int argc, char *argv[]) {
  // get file from command line
  Options opts(argc, argv);
  std::string parameters = opts.get_parameter_file();

  // Create a root
  pt::ptree root;

  // Load the ini file in this ptree
  pt::read_ini(parameters, root);

  double relerr_omega = root.get<double>("Friction.relerr_omega");
  // define needed quantities
  GreensTensor *greens_tensor = GreensTensorFactory::create(parameters);
  Polarizability *polarizabilty = PolarizabilityFactory::create(parameters);
  PowerSpectrumHarmOsc *powerspectrum =
      new PowerSpectrumHarmOsc(greens_tensor, polarizabilty);
  QuantumFriction *quant_friction = new QuantumFriction(
      greens_tensor, polarizabilty, powerspectrum, relerr_omega);

  // define looper
  Looper *looper = LooperFactory::create(parameters, quant_friction);

  // define output file
  std::ofstream file;
  file.open(opts.get_output_file());

  // define progressbar
  ProgressBar progbar(looper->get_steps_total(), 70);
  progbar.display();

  // calculate values
  double step, value;
  for (int i = 0; i < looper->get_steps_total(); i++) {
    step = looper->get_step(i);
    value = looper->calculate_value(i);

    ++progbar;
    progbar.display();
    file << step << "," << value << "\n";
  };

  // close file
  file.close();

  // close progress bar
  progbar.done();

  return 0;
};
