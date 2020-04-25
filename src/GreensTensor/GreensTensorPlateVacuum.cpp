// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "../ReflectionCoefficients/ReflectionCoefficientsFactory.h"
#include "GreensTensorPlate.h"
#include "GreensTensorPlateVacuum.h"

GreensTensorPlateVacuum::GreensTensorPlateVacuum(
    double v, double beta, double za,
    ReflectionCoefficients *reflection_coefficients, double delta_cut,
    vec::fixed<2> rel_err)
    : GreensTensorPlate(v, beta, za, reflection_coefficients, delta_cut,
                        rel_err) {
  this->vacuum_greens_tensor = new GreensTensorVacuum(v, beta, rel_err(0));
}

GreensTensorPlateVacuum::GreensTensorPlateVacuum(const std::string &input_file)
    : GreensTensorPlate(input_file) {
  // Create a root
  pt::ptree root;

  // Load the json file in this ptree
  pt::read_json(input_file, root);

  std::string addvacuum = root.get<std::string>("GreensTensor.addvacuum");
  assert(addvacuum == "true");

  this->vacuum_greens_tensor =
      new GreensTensorVacuum(v, beta, this->rel_err(0));
}

void GreensTensorPlateVacuum::integrate_k(cx_mat::fixed<3, 3> &GT,
                                          Options_GreensTensor opts) {

  cx_mat::fixed<3, 3> vac;

  GreensTensorPlate::integrate_k(GT, opts);

  opts.class_pt = vacuum_greens_tensor;
  vacuum_greens_tensor->integrate_k(vac, opts);

  GT += vac;
}

void GreensTensorPlateVacuum::calculate_tensor(cx_mat::fixed<3, 3> &GT,
                                               Options_GreensTensor opts) {

  cx_mat::fixed<3, 3> vac;

  GreensTensorPlate::calculate_tensor(GT, opts);

  opts.class_pt = vacuum_greens_tensor;
  vacuum_greens_tensor->calculate_tensor(vac, opts);

  GT += vac;
}

void GreensTensorPlateVacuum::print_info(std::ofstream &file) {
  file << "# GreensTensorPlateVacuum\n"
       << "# v = " << v << "\n"
       << "# beta = " << beta << "\n"
       << "# z_a = " << za << "\n"
       << "# delta_cut = " << delta_cut << "\n"
       << "# rel_err = " << rel_err(0) << " " << rel_err(1) << "\n";

  reflection_coefficients->print_info(file);
}
