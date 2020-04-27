// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "LooperV.h"
#include "../Friction/Friction.h"

LooperV::LooperV(double start, double end, int number_of_steps,
                 std::string scale)
    : Looper(start, end, number_of_steps, scale){};

LooperV::LooperV(std::string input_file) : Looper(input_file) {
  // Create a root
  pt::ptree root;

  // Load the json file in this ptree
  pt::read_json(input_file, root);

  // check if type is right
  std::string type = root.get<std::string>("Looper.type");
  assert(type == "v");
};

double LooperV::calculate_value(int step, void *quantity) {

  Friction *quantum_friction = static_cast<Friction *>(quantity);

  Options_Friction opts;
  opts.spectrum = NON_LTE_ONLY;
  // opts.full_spectrum = true;
  opts.class_pt = quantum_friction;

  // change v
  quantum_friction->get_greens_tensor()->set_v(this->steps[step]);

  return quantum_friction->calculate(opts);
};
