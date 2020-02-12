#include "GreensTensorPlateVacuum.h"
#include "../ReflectionCoefficients/ReflectionCoefficientsFactory.h"

GreensTensorPlateVacuum::GreensTensorPlateVacuum(
    double v, double za, double beta,
    ReflectionCoefficients *reflection_coefficients, double delta_cut,
    vec::fixed<2> rel_err)
    : GreensTensorPlate(v, za, beta, reflection_coefficients, delta_cut,
                        rel_err) {
  this->vacuum_greens_tensor = new GreensTensorVacuum(v, beta, rel_err(0));
};

GreensTensorPlateVacuum::GreensTensorPlateVacuum(std::string input_file)
    : GreensTensorPlate(input_file) {
  // Create a root
  pt::ptree root;

  // Load the ini file in this ptree
  pt::read_ini(input_file, root);

  std::string addvacuum = root.get<std::string>("GreensTensor.addvacuum");
  assert(addvacuum == "true");

  this->vacuum_greens_tensor =
      new GreensTensorVacuum(v, beta, this->rel_err(0));
};

void GreensTensorPlateVacuum::integrate_1d_k(cx_mat::fixed<3, 3> &GT,
                                             Options_GreensTensor opts) {

  cx_mat::fixed<3, 3> vac;

  GreensTensorPlate::integrate_1d_k(GT, opts);

  opts_vacuum = opts;
  opts.class_pt = vacuum_greens_tensor;
  vacuum_greens_tensor->integrate_1d_k(vac, opts);

  GT += vac;
};

void GreensTensorPlateVacuum::integrate_2d_k(cx_mat::fixed<3, 3> &GT,
                                             Options_GreensTensor opts) {

  cx_mat::fixed<3, 3> vac;

  GreensTensorPlate::integrate_2d_k(GT, opts);

  opts_vacuum = opts;
  opts.class_pt = vacuum_greens_tensor;
  vacuum_greens_tensor->integrate_2d_k(vac, opts);

  GT += vac;
};

void GreensTensorPlateVacuum::calculate_tensor(cx_mat::fixed<3, 3> &GT,
                                               Options_GreensTensor opts) {

  cx_mat::fixed<3, 3> vac;

  GreensTensorPlate::calculate_tensor(GT, opts);

  opts_vacuum = opts;
  opts.class_pt = vacuum_greens_tensor;
  vacuum_greens_tensor->calculate_tensor(vac, opts);

  GT += vac;
};
