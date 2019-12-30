// integration routine
#include "../Calculations/Integrations.h"

// ini parser
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
namespace pt = boost::property_tree;

#include "GreensTensorPlate.h"
#include "Permittivity/PermittivityFactory.h"


GreensTensorPlate::GreensTensorPlate(double v, double za, double beta, std::string input_file)
{
  // set velocity
  this->v  = v;
  // set distance from the surface
  this->za = za;
  // set inverse temperature
  this->beta = beta;

  // Create a root
  pt::ptree root;

  // Load the ini file in this ptree
  pt::read_ini(input_file, root);

  // read permittivity
  this->permittivity = PermittivityFactory::create(input_file);
};


std::complex<double> GreensTensorPlate::get_epsilon(double omega)
{
    return this->permittivity->epsilon(omega);
};

void GreensTensorPlate::calculate_tensor(cx_mat::fixed<3,3>& GT, vec::fixed<2> kvec, double omega)
{
  // calculating the solely the imaginary part of the free Green tensor
  double pre, kx, ky, k2, omega2;
  kx = kvec(0);
  ky = kvec(1);
  k2 = kx*kx + ky*ky;
  omega2 = omega*omega;
  GT.zeros();

//  GT(0,0) = pre*(omega2 - kx*kx);
//  GT(1,1) = pre*(omega2 - ky*ky);
//  GT(2,2) = pre*k2;
};

void GreensTensorPlate::integrate_k_2d(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts)
{
};


void GreensTensorPlate::integrate_k_1d(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts)
{
};

double GreensTensorPlate::integrand_k_1d(double kv, void *opts)
{
double result;

return result;
};
