#ifndef LOOPER_H
#define LOOPER_H

#include "../Friction/Friction.h"
#include <string>

class Looper {
protected:

  std::string scale; // scale type

  double start;              // starting value
  double end;                // end value
  int number_of_steps;       // number of steps
  std::vector<double> steps; // array containing the steps

  void calculate_steps();

public:
  // constructors
  Looper(double start, double end, int number_of_steps, std::string scale);
  Looper(std::string input_file);

  // calculate the the value of quantum friction
  virtual double calculate_value(int step, void *quantity) = 0;

  // getter functions
  int get_steps_total() { return this->number_of_steps; };
  double get_step(int i) { return this->steps[i]; };
};

#endif // LOOPER_H
