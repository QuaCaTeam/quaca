#include "LooperV.h"

LooperV::LooperV(double start, double end, int number_of_steps,
                 std::string scale, std::string type,
                 QuantumFriction *quantum_friction)
    : Looper(start, end, number_of_steps, scale, type, quantum_friction){};

LooperV::LooperV(std::string input_file, QuantumFriction *quantum_friction)
    : Looper(input_file, quantum_friction){};

double LooperV::calculate_value(int step) {
  Options_Friction opts;
  opts.class_pt = this->quantum_friction;

  // change v
  this->quantum_friction->greens_tensor->set_v(this->steps[step]);

  return this->quantum_friction->calculate(opts, 0., 1E3, 1e-5, 0);
};
