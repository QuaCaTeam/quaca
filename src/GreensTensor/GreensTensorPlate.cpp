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
  // read parameters
  this->delta_cut = root.get<double>("Plate_Setup.delta_cut");

};


std::complex<double> GreensTensorPlate::get_epsilon(double omega)
{
    return this->permittivity->epsilon(omega);
};

void GreensTensorPlate::calculate_tensor(cx_mat::fixed<3,3>& GT, vec::fixed<2> kvec, double omega)
{
  double pre, kx, ky, k_quad, omega_quad;
  // imaginary unit
  std::complex<double> I(0.0, 1.0);

  std::complex<double> eps, kappa, kappa_epsilon;
  std::complex<double> r_p, r_s;
  std::complex<double> prefactor_p, prefactor_s;


  kx = kvec(0);
  ky = kvec(1);
  k_quad = kx*kx + ky*ky;
  omega_quad = omega*omega;
  eps = permittivity->epsilon(omega);

  // kapppa as well as kappa_epsilon are defined to have either a purely
  // positive real part or purely negatively imaginary part

  kappa = sqrt(k_quad - omega_quad);
  kappa = ( std::abs(kappa.real()) ,
           -std::abs(kappa.imag()) );

  kappa_epsilon = sqrt(k_quad - eps*omega_quad);
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
  prefactor_s =  prefactor_p * r_s * omega_quad/(k_quad-omega_quad);
  prefactor_p = prefactor_p * r_p;

  // In the following, we already omit odd orders of ky, as we use that the
  // investigated system is symmetric in y and thus those terms vanish

  GT.zeros();
  GT(0,0) = prefactor_p * kx * kx / k_quad + prefactor_s * kx * kx / k_quad;
  GT(1,1) = prefactor_p * ky * ky / k_quad + prefactor_s * kx * kx / k_quad;
  GT(2,2) = prefactor_p * k_quad / kappa;
  GT(2,0) = I * prefactor_p * kx / kappa;
  GT(0,2) = - GT(2,0);
};



void GreensTensorPlate::integrate_k_1d(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts)
{

  // imaginary unit
  std::complex<double> I(0.0, 1.0);

  GT.zeros();
  opts.indices = {0,0};
  GT(0,0) = cquad(&integrand_k_1d, &opts, 0, M_PI ,1E-4,0) / M_PI;
  opts.indices = {1,1};
  GT(1,1) = cquad(&integrand_k_1d, &opts, 0, M_PI ,1E-4,0) / M_PI;
  opts.indices = {2,2};
  GT(2,2) = cquad(&integrand_k_1d, &opts, 0, M_PI ,1E-4,0) / M_PI;
  opts.indices = {2,0};
  GT(2,0) = I * cquad(&integrand_k_1d, &opts, 0, M_PI ,1E-4,0) / M_PI;
  GT(0,2) = - GT(2,0);
};

double GreensTensorPlate::integrand_k_1d(double phi, void *opts)
{
  //Units: c=1, 4 pi epsilon_0 = 1, hbar = 1
  Options_GreensTensor* opts_pt = static_cast<Options_GreensTensor*>(opts);
  GreensTensorPlate* pt = static_cast<GreensTensorPlate*>(opts_pt->class_pt);

  double omega = opts_pt->omega;

  double beta = pt->beta;
  double v = pt->v;
  double za = pt->za;
  double kappa_cut = pt->delta_cut / (2*za);
  double result;

  // Write the integration variable into the options struct
  opts_pt->kvec(0) = phi;

  // Calculate the integrand corresponding to the given options
  result = cquad(&integrand_k_2d, opts, - omega, kappa_cut ,1E-6,0);
  return result;
};

void GreensTensorPlate::integrate_k_2d(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts)
{
};
double GreensTensorPlate::integrand_k_2d(double kappa_double, void *opts)
{
  // This function calculates the integration with respect to kappa
  Options_GreensTensor* opts_pt = static_cast<Options_GreensTensor*>(opts);
  GreensTensorPlate* pt = static_cast<GreensTensorPlate*>(opts_pt->class_pt);

  // read omega and phi from the option struct
  double omega = opts_pt->omega;
  double cos_phi = cos(opts_pt->kvec(0));
  // read general input parameters
  double beta = pt->beta;
  double v = pt->v;
  double za = pt->za;

  // Before we can calculate the real or imaginary part of the chosen matrix
  // element, we need to store the complex result in result_complex.
  std::complex<double> result_complex;

  double result, omega_pl, omega_pl_quad, k_quad, kappa_quad;
  double omega_quad = omega * omega;
  double cos_phi_quad, sin_phi_quad, v_quad, k;
  cos_phi_quad = cos_phi * cos_phi;
  sin_phi_quad = 1.0 - cos_phi_quad;
  v_quad = v*v;

  // imaginary unit
  std::complex<double> I(0.0, 1.0);

  // permittivity and propagation through vacuum (kappa) and material (kappa_epsilon)
  std::complex<double> eps, kappa_complex, kappa_epsilon;

  // reflection coefficients and pre-factors of the corresponding polarization
  std::complex<double> r_p, r_s;
  std::complex<double> prefactor_p, prefactor_s;

  // Transfer kappa to the correct complex value
  if (kappa_double < 0.0) {
    kappa_complex = (0.0 , kappa_double);
    kappa_quad = - kappa_double * kappa_double ;
  }
  else
  {
    kappa_complex = (kappa_double , 0.0);
    kappa_quad =  kappa_double * kappa_double ;
  }

  k = (sqrt(kappa_quad*(1.0-v_quad*cos_phi_quad) + omega_quad) + v*omega*cos_phi)
  /(1 - v_quad*cos_phi_quad);

  k_quad = k*k;
  omega_pl = (omega+k * cos_phi * v);
  omega_pl_quad = omega_pl*omega_pl;

  eps = pt->get_epsilon(omega_pl);

  // kapppa as well as kappa_epsilon are defined to have either a purely
  // positive real part or purely negatively imaginary part
  kappa_epsilon = sqrt(k_quad - eps*omega_pl_quad);
  kappa_epsilon = ( std::abs(kappa_epsilon.real()) ,-std::abs(kappa_epsilon.imag()) );

  // Defining the reflection coefficients in transverse magnetice polarization (p)
  // and in transverse electric polarization (s)

  r_p = (eps*kappa_complex - kappa_epsilon)/
  (eps*kappa_complex + kappa_epsilon);

  r_s = (kappa_complex - kappa_epsilon)/
  (kappa_complex + kappa_epsilon);

  // For an better overview and a efficient calculation, we collect the pre-factors
  // of the p and s polarization separately

  prefactor_p =  exp(-2 * za * kappa_complex) /  (1. - cos_phi*v*omega_pl/k) ;
  prefactor_s =  prefactor_p * r_s * omega_pl_quad;
  prefactor_p =  prefactor_p * r_p * kappa_quad;

  // Calculate the G_xx element
  if(opts_pt->indices(0) == 0 && opts_pt->indices(1) == 0)
  {
    result_complex = prefactor_p * cos_phi_quad
    + prefactor_s * sin_phi_quad;
  }
  // Calculate the G_yy element
  else if(opts_pt->indices(0) == 1 && opts_pt->indices(1) == 1)
  {
    result_complex = prefactor_p * sin_phi_quad
    + prefactor_s * cos_phi_quad;
  }
  // Calculate the G_zz element
  else if(opts_pt->indices(0) == 2 && opts_pt->indices(1) == 2)
  {
    result_complex = prefactor_p *  k_quad / kappa_quad;
  }
  // Calculate the G_zx element
  else if(opts_pt->indices(0) == 2 && opts_pt->indices(1) == 0)
  {
    result_complex = I * prefactor_p * k * cos_phi / kappa_complex;
  }
  // Calculate the G_xz element
  else if(opts_pt->indices(0) == 0 && opts_pt->indices(1) == 2)
  {
    result_complex = - I * prefactor_p * k * cos_phi / kappa_complex;
  }
  else
  {
    result_complex = (0,0);
  }

  // Add weighting function if demanded
  if (opts_pt->fancy_I_kv)
  {
    result *= k * cos_phi;
  }
  else if(opts_pt->fancy_I_temp)
  {
    result /= (1.0-exp(-beta*omega_pl));
  }
  else if (opts_pt->fancy_I_kv_temp)
  {
    result *= k * cos_phi / (1.0-exp(-beta*omega_pl));
  }

  // Calculate fancy real part of the given matrix element
  if(opts_pt->fancy_R)
  {
    if (opts_pt->indices(0) == 2 && opts_pt->indices(1) == 0 ||
        opts_pt->indices(0) == 0 && opts_pt->indices(1) == 2 )
    {
      // Mind the missing leading I! This must be added after the double integration!
      result = result_complex.imag();
    }
    else
    {
      result = result_complex.real();
    }
  }
  // Calculate fancy imaginary part of the given matrix element
  else
  {
    if (opts_pt->indices(0) == 2 && opts_pt->indices(1) == 0 ||
        opts_pt->indices(0) == 0 && opts_pt->indices(1) == 2 )
    {
      // Mind the missing leading I! This must be added after the double integration!
      result = -result_complex.real();
    }
    else
    {
      result = result_complex.imag();
    }
  }
  return result;
};
