#ifndef LOOPERV_H
#define LOOPERV_H

#include "../Friction/Friction.h"
#include "Looper.h"
#include <string>

class LooperV : public Looper {
public:
  // constructors
  LooperV(double start, double end, int number_of_steps, std::string scale,
          Friction *quantum_friction);
  LooperV(std::string input_file, Friction *quantum_friction);

  // calculate the the value of quantum friction
  double calculate_value(int step);
};

#endif // LOOPERV_H
