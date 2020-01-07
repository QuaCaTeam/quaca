#include <assert.h>

// ini parser
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
namespace pt = boost::property_tree;

#include "GreensTensorVacuum.h"
#include "../Calculations/Integrations.h"

GreensTensorVacuum::GreensTensorVacuum(double v, double beta):GreensTensor(v,beta), za(za)
{
};

GreensTensorVacuum::GreensTensorVacuum(std::string input_file):GreensTensor(input_file)
{
};

void GreensTensorVacuum::calculate_tensor(cx_mat::fixed<3,3>& GT,Options_GreensTensor opts)
{
  // calculating the solely the imaginary part of the free Green tensor
  double pre, k_x, k_y, k_quad, omega, omega_quad;
  k_x = opts.kvec(0);
  k_y = opts.kvec(1);
  omega = opts.omega;
  k_quad = k_x*k_x + k_y*k_y;
  omega_quad = omega*omega;
  pre = 1.0/(2*M_PI*sqrt(omega_quad - k_quad));
  GT.zeros();
  GT(0,0) = pre*(omega_quad - k_x*k_x);
  GT(1,1) = pre*(omega_quad - k_y*k_y);
  GT(2,2) = pre*k_quad;
};

void GreensTensorVacuum::integrate_2d_k(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts)
{
  double k = opts.kvec(0);
  double omega = opts.omega;
  if(opts.fancy_I_kv)
  {
    // calculating the solely the imaginary part of the free Green tensor with Doppler shift in the frequency argument G(k, kx*v + omega)
    double xi_quad, omega_pl_quad, omega_quad;
    omega_quad = omega*omega;
    omega_pl_quad = (omega+k*this->v)*(omega+k*this->v);
    xi_quad = omega_pl_quad - k*k;
    GT.zeros();
    GT(0,0) = 0.5*xi_quad;
    GT(1,1) = 0.5*(omega_pl_quad - xi_quad*0.5);
    GT(2,2) = 0.5*(omega_pl_quad - xi_quad*0.5);
  }
  else
   {
     std::cerr << "Only the imaginary part of the Green tensor with Doppler shift in the frequency argument is implemented" << std::endl;
     exit(0);
   }

};


void GreensTensorVacuum::integrate_1d_k(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts)
{
  double omega = opts.omega;
  GT.zeros();

  if (opts.fancy_R)
  {
    std::cerr<< "Real part of the vacuum Green's is divergent and already included within the omega_a" << std::endl;
    exit(0);
  }

  opts.indices(0) = 0;
  opts.indices(1) = 0;
  GT(0,0) = cquad(&integrand_1d_k,&opts,-omega/(1.0+this->v),omega/(1.0-this->v),1E-7,0);
  opts.indices(0) = 1;
  opts.indices(1) = 1;
  GT(1,1) = cquad(&integrand_1d_k,&opts, -omega/(1.0+this->v),omega/(1.0-this->v),1E-7,0);
  GT(2,2) = GT(1,1);
};

double GreensTensorVacuum::integrand_1d_k(double kv, void *opts)
{
  //Units: c=1, 4 pi epsilon_0 = 1, hbar = 1
  Options_GreensTensor* opts_pt = static_cast<Options_GreensTensor*>(opts);
  GreensTensorVacuum* pt = static_cast<GreensTensorVacuum*>(opts_pt->class_pt);

  double omega = opts_pt->omega;
  double beta = pt->beta;
  double result, xi_quad, omega_pl, omega_pl_quad;

  omega_pl = (omega+kv*pt->v);
  omega_pl_quad = omega_pl*omega_pl;
  xi_quad = omega_pl_quad - kv*kv;

  if(opts_pt->indices(0) == 0 && opts_pt->indices(1) == 0)
  {
    result = 0.5*xi_quad;
  }
  else if(opts_pt->indices(0) == opts_pt->indices(1))
  {
    result = 0.5*(omega_pl_quad - xi_quad*0.5);
  }
  else
  {
    return 0;
  }
  if(opts_pt->fancy_I)
  {
    return result;
  }
  if (opts_pt->fancy_I_kv)
  {
    result *= kv;
  }
  else if(opts_pt->fancy_I_temp)
  {
    result /= (1.0-exp(-beta*omega_pl));
  }
  else if (opts_pt->fancy_I_kv_temp)
  {
    result *= kv/(1.0-exp(-beta*omega_pl));
  }
  return result;
};
