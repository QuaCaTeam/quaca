#ifndef GREENSTENSORPLATE_H
#define GREENSTENSORPLATE_H

#include "GreensTensor.h"

class GreensTensorPlate : public GreensTensor
{
public:

  // constructors
  GreensTensorPlate(double v, double za, double beta): GreensTensor(v, za, beta) {};
  GreensTensorPlate(std::string input_file): GreensTensor(input_file) {};

  // calculate full tensor
  void calculate_tensor(cx_mat::fixed<3,3>& GT, vec::fixed<2> kvec, double omega);

  // integrate over a two-dimensional k space
  void integrate_2d_k(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts);

  // integrate over a one-dimensional k space
  void integrate_1d_k(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts);

};


#endif // GREENSTENSORPLATE_H
