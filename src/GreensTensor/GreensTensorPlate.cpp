#include "GreensTensorPlate.h"

// ini parser
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

GreensTensorPlate::GreensTensorPlate(std::string input_file)
    : GreensTensor(input_file) {
  pt::ptree root;
  pt::read_ini(input_file, root);
  this->z_a = root.get<double>("GreensTensor.z_a");
};

GreensTensorPlate::GreensTensorPlate(double v, double z_a, double beta)
    : GreensTensor(v, beta), z_a(z_a) {
  assert(z_a > 0);
};

void GreensTensorPlate::calculate_tensor(cx_mat::fixed<3, 3> &GT,
                                         Options_GreensTensor opts) {
  double a = 0;
};

void GreensTensorPlate::integrate_2d_k(cx_mat::fixed<3, 3> &GT,
                                       Options_GreensTensor opts) {
  double a = 0;
};

void GreensTensorPlate::integrate_1d_k(cx_mat::fixed<3, 3> &GT,
                                       Options_GreensTensor opts) {
  double a = 0;
};
