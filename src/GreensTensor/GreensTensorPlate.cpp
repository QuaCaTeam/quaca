// ini parser
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
namespace pt = boost::property_tree;

#include "GreensTensorPlate.h"

void GreensTensorPlate::calculate_tensor(cx_mat::fixed<3,3>& GT, vec::fixed<2> kvec, double omega)
{
  double a = 0;
};

void GreensTensorPlate::integrate_2d_k(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts)
{
  double a = 0;
};


void GreensTensorPlate::integrate_1d_k(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts)
{
  double a = 0;
};
