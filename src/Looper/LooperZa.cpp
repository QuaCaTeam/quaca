// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <utility>
namespace pt = boost::property_tree;

#include "../GreensTensor/GreensTensorPlate.h"
#include "LooperZa.h"

LooperZa::LooperZa(double start, double end, int number_of_steps,
                   std::string scale)
    : Looper(start, end, number_of_steps, std::move(scale)) {}

LooperZa::LooperZa(const std::string &input_file) : Looper(input_file) {
  // Create a root
  pt::ptree root;

  // Load the json file in this ptree
  pt::read_json(input_file, root);

  // check if type is right
  std::string type = root.get<std::string>("Looper.type");
  assert(type == "za");
}

double LooperZa::calculate_value(int step, Friction *quantum_friction) {
  Options_Friction opts;
  opts.spectrum = NON_LTE_ONLY;
  opts.class_pt = quantum_friction;

  // change za
  auto *pt =
      dynamic_cast<GreensTensorPlate *>(quantum_friction->get_greens_tensor());

  if (pt == nullptr) {
    std::cerr
        << "You try to loop over z_a, but did not give a plate Green's tensor."
        << std::endl;
    exit(-1);
  }

  pt->set_za(this->steps[step]);

  return quantum_friction->calculate(opts);
}
