#ifndef LOOPERZA_H
#define LOOPERZA_H

#include "../Calculations/QuantumFriction.h"
#include "Looper.h"
#include <string>

class LooperZa : public Looper {
public:
  // constructors
  LooperZa(double start, double end, int number_of_steps, std::string scale,
           QuantumFriction *quantum_friction);
  LooperZa(std::string input_file, QuantumFriction *quantum_friction);

  // calculate the the value of quantum friction
  double calculate_value(int step);
};

#endif // LOOPERZA_H
