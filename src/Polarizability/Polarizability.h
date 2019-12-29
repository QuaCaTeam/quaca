#ifndef POLARIZABILITY_H
#define POLARIZABILITY_H

#include <complex>
#include <cmath>
#include <armadillo>

#include "../GreensTensor/GreensTensor.h"
#include "MemoryKernel/MemoryKernel.h"

using namespace arma;

//! a struct with the integration options
struct Options_Polarizabiliy;

//! An abstract polarizability class
/*!
* This is an abstract polarizability class.
* All polarizabilities should return a 3 by 3 matrix, given a real frequency as input.
* The only two child class of this will be 1) a polarizability where the particle
* has no interval bath 2) a polarizability where the particle has an internal bath
*/
class Polarizability
{
protected:

    double omega_a;              // resonance frequency
    double alpha_zero;           // vacuum polarizability
    GreensTensor *greens_tensor; // Green's tensor

public:

    // constructors
    Polarizability(double omega_a, double alpha_zero, GreensTensor *greens_tensor): omega_a(omega_a), alpha_zero(alpha_zero), greens_tensor(greens_tensor) {};
    Polarizability(std::string input_file);

    // calculate the polarizability tensor
    virtual void calculate(cx_mat::fixed<3,3>& alpha,double omega) =0;

    // getter functions
    double get_omega_a(){return this->omega_a;};
    double get_alpha_zero(){return this->alpha_zero;};

};

// A struct for integration options
struct Options_Polarizability
{
 // Different options for the integrand
  bool fancy_R = false;
  bool fancy_I = false;
  bool fancy_I_kv = false;
  bool fancy_I_temp = false;
  bool fancy_I_kv_temp = false;
  //Indices of the 3x3 polarizability tensor
  vec::fixed<2> indices = {-1,-1};
  //Value of omega for the integration of the k-Variables
  double omega = NAN;
  //Pointer to the Polarizability to be able to access the attributes of the class eventhough the integrand is static
  Polarizability* class_pt;
};

#endif //POLARIZABILITY_H
