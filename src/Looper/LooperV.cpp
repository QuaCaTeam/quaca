// ini parser
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "LooperV.h"

LooperV::LooperV(double start, double end, int number_of_steps,
                 std::string scale, QuantumFriction *quantum_friction)
    : Looper(start, end, number_of_steps, scale, quantum_friction){};

LooperV::LooperV(std::string input_file, QuantumFriction *quantum_friction)
    : Looper(input_file, quantum_friction) {
  // Create a root
  pt::ptree root;

  // Load the ini file in this ptree
  pt::read_ini(input_file, root);

  // check if type is right
  std::string type = root.get<std::string>("Looper.type");
  assert(type == "v");
};

double LooperV::calculate_value(int step) {
  Options_Friction opts;
  opts.non_LTE = true;
  // opts.full_spectrum = true;
  opts.class_pt = this->quantum_friction;

  // change v
  this->quantum_friction->greens_tensor->set_v(this->steps[step]);

  return this->quantum_friction->calculate(opts, 1e-1, 0);
};
