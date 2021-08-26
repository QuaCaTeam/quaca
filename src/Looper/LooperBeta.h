#ifndef LOOPERBETA_H
#define LOOPERBETA_H

#include "../Friction/Friction.h"
#include "Looper.h"
#include <string>

class LooperBeta : public Looper {
public:
  // constructors
  LooperBeta(double start, double end, int number_of_steps, const std::string &scale);
  LooperBeta(const std::string &input_file);

  // calculate the the value of quantum friction
  double calculate_value(int step, std::shared_ptr<Friction> quantum_friction) const override;

  // print info
  void print_info(std::ostream &stream) const override;
};

#endif // LOOPERBETA_H
