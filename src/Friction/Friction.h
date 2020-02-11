#ifndef QUANTUMFRICTION_H
#define QUANTUMFRICTION_H

#include "../GreensTensor/GreensTensor.h"
#include "../Polarizability/Polarizability.h"
#include "../PowerSpectrum/PowerSpectrum.h"
#include <string>

/*!
 * This is a class computing the quantum friction force for a given Green's
 * tensor and polarizibility with a certain set of parameters
 */

// declaration of struct containing all options for the calculation of the
// quantum friction force, definition can be found below
struct Options_Friction;

class Friction {
protected:
  double relerr_omega;

  GreensTensor *greens_tensor;
  Polarizability *polarizability;
  PowerSpectrum *powerspectrum;

public:
  Friction(std::string input_file);
  Friction(GreensTensor *greens_tensor, Polarizability *polarizability,
                  PowerSpectrum *powerspectrum, double relerr_omega);
  double calculate(Options_Friction opts);
  static double friction_integrand(double omega, void *opts);

  // getter functions
  GreensTensor *get_greens_tensor() { return greens_tensor; };
  Polarizability *get_polarizability() { return polarizability; };
  PowerSpectrum *get_powerspectrum() { return powerspectrum; };
};

struct Options_Friction {
  Friction *class_pt;
  bool full_spectrum = false;
  bool non_LTE = false;
};

#endif // QUANTUMFRICTION_H
