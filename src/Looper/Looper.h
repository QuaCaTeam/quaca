#ifndef LOOPER_H
#define LOOPER_H

#include "../Calculations/QuantumFriction.h"
#include <string>

class Looper {
protected:
  std::string scale; // scale type
  std::string type;  // type of looper, that is which variable gets iterated

  double start;              // starting value
  double end;                // end value
  int number_of_steps;       // number of steps
  std::vector<double> steps; // array containing the steps

  void calculate_steps();

public:
  QuantumFriction *quantum_friction;

  // constructors
  Looper(double start, double end, int number_of_steps, std::string scale,
         std::string type, QuantumFriction *quantum_friction);
  Looper(std::string input_file, QuantumFriction *quantum_friction);

  // calculate the the value of quantum friction
  virtual double calculate_value(int step) = 0;

  // getter functions
  int get_steps_total() { return this->number_of_steps; };
  double get_step(int i) { return this->steps[i]; };
};

#endif // LOOPER_H
