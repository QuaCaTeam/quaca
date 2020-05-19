#ifndef LOOPERV_H
#define LOOPERV_H

#include "../Friction/Friction.h"
#include "Looper.h"
#include <string>

class LooperV : public Looper {
public:
  // constructors
  LooperV(double start, double end, int number_of_steps, const std::string &scale);
  LooperV(const std::string &input_file);

  // calculate the the value of quantum friction
  double calculate_value(int step, std::shared_ptr<Friction> quantum_friction) const override;

  // print info
  void print_info(std::ostream &stream) const override;
};

#endif // LOOPERV_H
