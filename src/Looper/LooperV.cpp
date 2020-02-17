// ini parser
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "LooperV.h"

LooperV::LooperV(double start, double end, int number_of_steps,
                 std::string scale)
      : Looper(start, end, number_of_steps, scale){};

LooperV::LooperV(std::string input_file): Looper(input_file){
  // Create a root
  pt::ptree root;

  // Load the ini file in this ptree
  pt::read_ini(input_file, root);

  // check if type is right
  std::string type = root.get<std::string>("Looper.type");
  assert(type == "v");
};

double LooperV::calculate_value(int step, Friction* quantum_friction) {
  Options_Friction opts;
  opts.non_LTE = true;
  // opts.full_spectrum = true;
  opts.class_pt = quantum_friction;

  // change v
  quantum_friction->get_greens_tensor()->set_v(this->steps[step]);

  return quantum_friction->calculate(opts);
};
