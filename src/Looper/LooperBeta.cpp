// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "LooperBeta.h"

LooperBeta::LooperBeta(double start, double end, int number_of_steps,
                 const std::string &scale)
    : Looper(start, end, number_of_steps, scale) {}

LooperBeta::LooperBeta(const std::string &input_file) : Looper(input_file) {
  // Create a root
  pt::ptree root;

  // Load the json file in this ptree
  pt::read_json(input_file, root);

  // check if type is right
  std::string type = root.get<std::string>("Looper.type");
  assert(type == "beta");
}

double
LooperBeta::calculate_value(int step,
                         std::shared_ptr<Friction> quantum_friction) const {
  // change beta
  // TODO: set temperature everywhere, where it is needed
  quantum_friction->get_greens_tensor()->set_beta(this->steps[step]);

  return quantum_friction->calculate(NON_LTE_ONLY);
}

void LooperBeta::print_info(std::ostream &stream) const {
  stream << "# LooperBeta\n#\n"
         << "# start = " << start << "\n"
         << "# end = " << end << "\n"
         << "# number_of_steps = " << number_of_steps << "\n"
         << "# scale = " << scale << "\n";
}
