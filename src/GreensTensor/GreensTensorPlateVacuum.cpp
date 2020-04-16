// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <utility>
namespace pt = boost::property_tree;

#include "../ReflectionCoefficients/ReflectionCoefficientsFactory.h"
#include "GreensTensorPlate.h"
#include "GreensTensorPlateVacuum.h"

GreensTensorPlateVacuum::GreensTensorPlateVacuum(
    double v, double beta, double za,
    std::shared_ptr<ReflectionCoefficients> reflection_coefficients,
    double delta_cut, vec::fixed<2> rel_err)
    : GreensTensorPlate(v, beta, za, std::move(reflection_coefficients),
                        delta_cut, rel_err) {
  this->vacuum_greens_tensor =
      std::make_shared<GreensTensorVacuum>(v, beta, rel_err(0));
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
      std::make_shared<GreensTensorVacuum>(v, beta, this->rel_err(0));
}

void GreensTensorPlateVacuum::integrate_k(cx_mat::fixed<3, 3> &GT,
                                          Options_GreensTensor opts) {

  cx_mat::fixed<3, 3> vac;

  GreensTensorPlate::integrate_k(GT, opts);

  opts_vacuum = opts;
  opts.class_pt = vacuum_greens_tensor;
  vacuum_greens_tensor->integrate_k(vac, opts);

  GT += vac;
}

void GreensTensorPlateVacuum::calculate_tensor(cx_mat::fixed<3, 3> &GT,
                                               Options_GreensTensor opts) {

  cx_mat::fixed<3, 3> vac;

  GreensTensorPlate::calculate_tensor(GT, opts);

  opts_vacuum = opts;
  opts.class_pt = vacuum_greens_tensor;
  vacuum_greens_tensor->calculate_tensor(vac, opts);

  GT += vac;
}
