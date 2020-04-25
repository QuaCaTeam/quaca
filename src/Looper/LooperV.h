#ifndef LOOPERV_H
#define LOOPERV_H

#include "../Friction/Friction.h"
#include "Looper.h"
#include <string>

class LooperV : public Looper {
public:
  // constructors
  LooperV(double start, double end, int number_of_steps, std::string scale);
  LooperV(const std::string& input_file);

  // calculate the the value of quantum friction
  double calculate_value(int step, Friction* quantum_friction);
};

#endif // LOOPERV_H
