#ifndef GREENSTENSORPLATE_H
#define GREENSTENSORPLATE_H

#include "GreensTensor.h"

class GreensTensorPlate : public GreensTensor {
private:
  double z_a; // distance from plate

public:
  // constructors
  GreensTensorPlate(std::string input_file);
  GreensTensorPlate(double v, double z_a, double beta);

  // calculate full tensor
  void calculate_tensor(cx_mat::fixed<3, 3> &GT, Options_GreensTensor opts);

  // integrate over a two-dimensional k space
  void integrate_2d_k(cx_mat::fixed<3, 3> &GT, Options_GreensTensor opts);

  // integrate over a one-dimensional k space
  void integrate_1d_k(cx_mat::fixed<3, 3> &GT, Options_GreensTensor opts);
};

#endif // GREENSTENSORPLATE_H
