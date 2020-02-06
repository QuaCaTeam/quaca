// ini parser
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "../GreensTensor/GreensTensorPlate.h"
#include "LooperZa.h"

LooperZa::LooperZa(double start, double end, int number_of_steps,
                   std::string scale, QuantumFriction *quantum_friction)
    : Looper(start, end, number_of_steps, scale, quantum_friction){};

LooperZa::LooperZa(std::string input_file, QuantumFriction *quantum_friction)
    : Looper(input_file, quantum_friction) {
  // Create a root
  pt::ptree root;

  // Load the ini file in this ptree
  pt::read_ini(input_file, root);

  // check if type is right
  std::string type = root.get<std::string>("Looper.type");
  assert(type == "za");
};

double LooperZa::calculate_value(int step) {
  Options_Friction opts;
  opts.non_LTE = true;
  opts.class_pt = this->quantum_friction;

  // change za
  GreensTensorPlate *pt =
      dynamic_cast<GreensTensorPlate *>(this->quantum_friction->greens_tensor);

  if (pt == NULL) {
    std::cerr
        << "You try to loop over z_a, but did not give a plate Green's tensor."
        << std::endl;
    exit(-1);
  }

  pt->set_z_a(this->steps[step]);

  return this->quantum_friction->calculate(opts, 1e-2, 0);
};
