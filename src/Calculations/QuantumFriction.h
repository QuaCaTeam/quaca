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
  double relerr_omega;

  GreensTensor *greens_tensor;
  Polarizability *polarizability;
  PowerSpectrum *powerspectrum;

  //Constructor with ini-file
  QuantumFriction(std::string input_file);
  //Constructor with initialization list
  QuantumFriction(GreensTensor *greens_tensor, Polarizability *polarizability,
                  PowerSpectrum *powerspectrum, double relerr_omega);
  double calculate(Options_Friction opts);
  static double friction_integrand(double omega, void *opts);
};

//Integration options for the quantum friction calculation
struct Options_Friction {
    //Pointer storing the address of the class itself (similar to 'this')
  QuantumFriction *class_pt;
  //Compute the full power spectrum, using eq. (4.3) in Marty PhD thesis
  bool full_spectrum = false;
  //Use eq. (4.5) in Marty's PhD thesis
  bool non_LTE = false;
};

#endif // QUANTUMFRICTION_H
