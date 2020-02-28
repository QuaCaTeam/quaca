#ifndef GREENSTENSORPLATEMAGNETIC_H
#define GREENSTENSORPLATEMAGNETIC_H

#include <armadillo>
#include <assert.h>
#include "GreensTensorPlate.h"

struct Options_GreensTensorMagnetic: Options_GreensTensor
{
    Tensor_Options EB = IGNORE;
    Tensor_Options BE = IGNORE;
    Tensor_Options BB = IGNORE;
};

class GreensTensorPlateMagnetic: public GreensTensorPlate{
public:
    GreensTensorPlateMagnetic(std::string input_file): GreensTensorPlate(input_file) {};
    GreensTensorPlateMagnetic(double v, double za, double beta,
                               ReflectionCoefficients *reflection_coefficients,
                               double delta_cut, vec::fixed<2> rel_err ): GreensTensorPlate(v,za,beta,
                                                                       reflection_coefficients, delta_cut, rel_err) {};
   // calculate the tensor in frequency and momentum space
   void calculate_tensor(cx_mat::fixed<3, 3> &GT, Options_GreensTensor opts);

    // integrate over a two-dimensional k space
    //void integrate_k(cx_mat::fixed<3, 3> &GT, Options_GreensTensor opts);

    // integrands
    //static double integrand_1d_k(double kx, void *opts);
    static double integrand_2d_k_magnetic(double ky, void *opts);
};

#endif
