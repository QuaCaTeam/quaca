#ifndef LOOPER_H
#define LOOPER_H

#include "../Friction/Friction.h"
#include <string>

class Looper {
protected:
  double start;              // starting value
  double end;                // end value
  int number_of_steps;       // number of steps
  std::string scale;         // scale type
  std::vector<double> steps; // array containing the steps

  void calculate_steps();

public:
  // constructors
  Looper(double start, double end, int number_of_steps,
         const std::string &scale);
  explicit Looper(const std::string &input_file);

  // calculate the the value of quantum friction
  virtual double
  calculate_value(int step,
                  std::shared_ptr<Friction> quantum_friction) const = 0;

  // getter functions
  int get_steps_total() const { return this->number_of_steps; };
  double get_step(int i) const { return this->steps[i]; };

  // print info
  virtual void print_info(std::ostream &stream) const = 0;
};

#endif // LOOPER_H
