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

class QuantumFriction {
public:
  GreensTensor *greens_tensor;
  Polarizability *polarizability;
  PowerSpectrum *powerspectrum;

  QuantumFriction(std::string input_file);
  QuantumFriction(GreensTensor *greens_tensor, Polarizability *polarizability,
                  PowerSpectrum *powerspectrum);
  double calculate(Options_Friction opts, double relerr, double epsabs);
  static double friction_integrand(double omega, void *opts);
};

struct Options_Friction {
  QuantumFriction *class_pt;
  bool full_spectrum = false;
  bool non_LTE = false;
};

#endif
