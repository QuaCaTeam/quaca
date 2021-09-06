// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "LooperV.h"

LooperV::LooperV(double start, double end, int number_of_steps,
                 const std::string &scale)
    : Looper(start, end, number_of_steps, scale) {}

LooperV::LooperV(const std::string &input_file) : Looper(input_file) {
  // Create a root
  pt::ptree root;

  // Load the json file in this ptree
  pt::read_json(input_file, root);

  // check if type is right
  std::string type = root.get<std::string>("Looper.type");
  assert(type == "v");
}

double
LooperV::calculate_value(int step,
                         std::shared_ptr<Friction> quantum_friction) const {
  // change v
  quantum_friction->get_greens_tensor()->set_v(this->steps[step]);

  return quantum_friction->calculate(NON_LTE_ONLY);
}

void LooperV::print_info(std::ostream &stream) const {
  stream << "# LooperV\n#\n"
         << "# start = " << start << "\n"
         << "# end = " << end << "\n"
         << "# number_of_steps = " << number_of_steps << "\n"
         << "# scale = " << scale << "\n";
}
