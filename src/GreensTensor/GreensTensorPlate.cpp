// ini parser
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
namespace pt = boost::property_tree;

#include "GreensTensorPlate.h"

GreensTensorPlate::GreensTensorPlate(double v, double za, double beta)
{
  // set velocity
  this->v = v;
  // set za
  this->za = za;
  // set inverse temperature
  this->beta = beta;
};

GreensTensorPlate::GreensTensorPlate(std::string input_file)
{
  // Create a root
  pt::ptree root;

  // Load the ini file in this ptree
  pt::read_ini(input_file, root);

  // read parameters
  this->v = root.get<double>("GreensTensor.v");
  this->za = root.get<double>("GreensTensor.za");
  this->beta = root.get<double>("GreensTensor.beta");
};


void GreensTensorPlate::calculate_tensor(cx_mat::fixed<3,3>& GT, vec::fixed<2> kvec, double omega)
{
  double a = 0;
};

void GreensTensorPlate::integrate_k_2d(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts)
{
  double a = 0;
};


void GreensTensorPlate::integrate_k_1d(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts)
{
  double a = 0;
};
