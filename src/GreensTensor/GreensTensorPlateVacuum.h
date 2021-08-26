#ifndef GREENSTENSORPLATEVACUUM_H
#define GREENSTENSORPLATEVACUUM_H

#include <memory>
#include "GreensTensorPlate.h"
#include "GreensTensorVacuum.h"

class GreensTensorPlateVacuum : public GreensTensorPlate {
private:
  std::shared_ptr<GreensTensorVacuum> vacuum_greens_tensor;

public:
  // constructors
  GreensTensorPlateVacuum(double v, double beta, double za,
                          std::shared_ptr<ReflectionCoefficients> reflection_coefficients,
                          double delta_cut, vec::fixed<2> rel_err);
  explicit GreensTensorPlateVacuum(const std::string &input_file);

  // calculate the tensor in frequency and momentum space
  void calculate_tensor(double omega, vec::fixed<2> k,
                        cx_mat::fixed<3, 3> &GT) const override;

  // integrate over a two-dimensional k space
  void integrate_k(double omega, cx_mat::fixed<3, 3> &GT,
                   Tensor_Options fancy_complex,
                   Weight_Options weight_function) const override;

  // getters
  std::shared_ptr<GreensTensorVacuum> get_vacuums_greens_tensor() {
    return vacuum_greens_tensor;
  };

  // setters
  void set_v(double v) override {
    this->v = v;
    this->vacuum_greens_tensor->set_v(v);
  };
  void set_beta(double beta) override {
    this->beta = beta;
    this->vacuum_greens_tensor->set_beta(beta);
  };

  // print info
  void print_info(std::ostream &stream) const override;
};

#endif // GREENSTENSORPLATEVACUUM_H
