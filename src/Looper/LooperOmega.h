#ifndef LOOPEROMEGA_H
#define LOOPEROMEGA_H

#include "../Polarizability/Polarizability.h"
#include "Looper.h"
#include <string>

class LooperOmega : public Looper {
public:
  // constructors
  LooperOmega(double start, double end, int number_of_steps, std::string scale);
  LooperOmega(std::string input_file);

  // calculate the the value of quantum friction
  double calculate_value(int step, void* quantity);
};

#endif // LOOPEROMEGA_H
