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
  // imaginary unit
  std::complex<double> I(0.0, 1.0);

  std::complex<double> eps, kappa, kappa_epsilon;
  std::complex<double> r_p, r_s;
  std::complex<double> prefactor_p, prefactor_s;


  kx = kvec(0);
  ky = kvec(1);
  k2 = kx*kx + ky*ky;
  omega2 = omega*omega;
  eps = permittivity->epsilon(omega);

  // kapppa as well as kappa_epsilon are defined to have either a purely
  // positive real part or purely negatively imaginary part

  kappa = sqrt(k2 - omega2);
  kappa = ( std::abs(kappa.real()) ,
           -std::abs(kappa.imag()) );

  kappa_epsilon = sqrt(k2 - eps*omega2);
  kappa_epsilon = ( std::abs(kappa_epsilon.real()) ,
                   -std::abs(kappa_epsilon.imag()) );

  // Defining the reflection coefficients in transverse magnetice polarization (p)
  // and in transverse electric polarization (s)

  r_p = (eps*kappa - kappa_epsilon)/
        (eps*kappa + kappa_epsilon);

  r_s = (kappa - kappa_epsilon)/
        (kappa + kappa_epsilon);

  // For an better overview and a efficient calculation, we collect the pre-factors
  // of the p and s polarization separately

  prefactor_p = 2 * M_PI * kappa * std::exp(-2 * za * kappa);
  prefactor_s =  prefactor_p * r_s * omega2/(k2-omega2);
  prefactor_p = prefactor_p * r_p;

  // In the following, we already omit odd orders of ky, as we use that the
  // investigated system is symmetric in y and thus those terms vanish

  GT.zeros();
  GT(0,0) = prefactor_p * kx * kx / k2 + prefactor_s * ky * ky / k2;
  GT(1,1) = prefactor_p * ky * ky / k2 + prefactor_s * kx * kx / k2;
  GT(2,2) = prefactor_p * k2 / (k2 - omega2);
  GT(2,0) = I * prefactor_p * kx / (k2 - omega2);
  GT(0,2) = - GT(2,0);
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
