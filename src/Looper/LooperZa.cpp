#include "LooperZa.h"
#include "../GreensTensor/GreensTensorPlate.h"

LooperZa::LooperZa(double start, double end, int number_of_steps,
                   std::string scale, std::string type,
                   QuantumFriction *quantum_friction)
    : Looper(start, end, number_of_steps, scale, type, quantum_friction){};

LooperZa::LooperZa(std::string input_file, QuantumFriction *quantum_friction)
    : Looper(input_file, quantum_friction){};

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
