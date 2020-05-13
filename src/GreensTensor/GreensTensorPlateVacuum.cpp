// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <utility>
namespace pt = boost::property_tree;

#include "../ReflectionCoefficients/ReflectionCoefficientsFactory.h"
#include "GreensTensorPlate.h"
#include "GreensTensorVacuum.h"
#include "GreensTensorPlateVacuum.h"

GreensTensorPlateVacuum::GreensTensorPlateVacuum(
    double v, double beta, double za,
    std::shared_ptr<ReflectionCoefficients> reflection_coefficients,
    double delta_cut, vec::fixed<2> rel_err)
    : GreensTensorPlate(v, beta, za, std::move(reflection_coefficients), delta_cut,
                        rel_err) {
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

void GreensTensorPlateVacuum::integrate_k(
    double omega, cx_mat::fixed<3, 3> &GT, 
    Tensor_Options fancy_complex, Weight_Options weight_function) const {

  //compute the contributions from the planar surface
  GreensTensorPlate::integrate_k(omega, GT, fancy_complex, weight_function);

  //compute the contributions from the vacuum
  cx_mat::fixed<3, 3> vac;
  vacuum_greens_tensor->integrate_k(omega, vac, fancy_complex, weight_function);

  GT += vac;
}

void GreensTensorPlateVacuum::calculate_tensor(double omega, vec::fixed<2> k,
                                               cx_mat::fixed<3, 3> &GT) const {

  //compute the contributions from the planar surface
  GreensTensorPlate::calculate_tensor(omega, k, GT);

  //compute the contributions from the vacuum
  cx_mat::fixed<3, 3> vac;
  vacuum_greens_tensor->calculate_tensor(omega, k, vac);

  GT += vac;
}
