#include "Quaca.h"
#include <iostream>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

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
  QuantumFriction *quant_friction =
      new QuantumFriction(greens_tensor, polarizabilty, powerspectrum, relerr_omega);

  // define output file
  std::ofstream file;
  file.open(opts.get_output_file());

  /*
   * calculate quantum friction for wanted parameters
   */

  // close file
  file.close();

  // delete dynamically allocated memory
  delete[] greens_tensor;
  delete[] polarizabilty;
  delete[] powerspectrum;
  delete[] quant_friction;

  return 0;
};
