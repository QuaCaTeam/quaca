#ifndef GREENSTENSORPLATEVACUUM_H
#define GREENSTENSORPLATEVACUUM_H

#include "GreensTensorPlate.h"
#include "GreensTensorVacuum.h"

class GreensTensorPlateVacuum : public GreensTensorPlate {
private:
  GreensTensorVacuum *vacuum_greens_tensor;
  Options_GreensTensor opts_vacuum;

public:
  // constructors
  GreensTensorPlateVacuum(double v, double za, double beta,
                          ReflectionCoefficients *reflection_coefficients,
                          double delta_cut, vec::fixed<2> rel_err);
  GreensTensorPlateVacuum(std::string input_file);

  // calculate the tensor in frequency and momentum space
  void calculate_tensor(cx_mat::fixed<3, 3> &GT, Options_GreensTensor opts);

  // integrate over a two-dimensional k space
  void integrate_k(cx_mat::fixed<3, 3> &GT, Options_GreensTensor opts);

  // getters
  GreensTensorVacuum *get_vacuums_greens_tensor() {
    return vacuum_greens_tensor;
  };

  // setters
  void set_v(double v) {
    this->v = v;
    this->vacuum_greens_tensor->set_v(v);
  };
};

#endif // GREENSTENSORPLATEVACUUM_H
