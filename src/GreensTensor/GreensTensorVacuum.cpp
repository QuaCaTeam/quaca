#include "GreensTensorVacuum.h"
#include "../Calculations/Integrations.h"



GreensTensorVacuum::GreensTensorVacuum(double v, double beta)
{
  // set velocity
  this->v = v;
  // set inverse temperature
  this->beta = beta;

};
void GreensTensorVacuum::calculate_tensor(cx_mat::fixed<3,3>& GT, vec::fixed<2> kvec, double omega)
{
  // calculating the solely the imaginary part of the free Green tensor
  double pre, kx, ky, k2, omega2;
  kx = kvec(0);
  ky = kvec(1);
  k2 = kx*kx + ky*ky;
  omega2 = omega*omega;
  pre = 1.0/(2*M_PI*sqrt(omega2 - k2));
  GT.zeros();
  GT(0,0) = pre*(omega2 - kx*kx);
  GT(1,1) = pre*(omega2 - ky*ky);
  GT(2,2) = pre*k2;
};
void GreensTensorVacuum::integrate_k_2d(cx_mat::fixed<3,3>& GT, double omega, double kv, Options opts)
{
  if(opts.fancy_I_kv)
  {
    // calculating the solely the imaginary part of the free Green tensor with Doppler shift in the frequency argument G(k, kx*v + omega)
    double pre, xi2, omegapl2;
    omegapl2 = (omega+kv*this->v)*(omega+kv*this->v);
    xi2 = omegapl2 - kv*kv;
    pre = 0.5;
    GT.zeros();
    GT(0,0) = pre*xi2;
    GT(1,1) = pre*(omegapl2 - xi2*0.5);
    GT(2,2) = pre*(omegapl2 - xi2*0.5);
  }
  else
   {
     std::cerr << "Only the imaginary part of the Green tensor with Doppler shift in the frequency argument is implemented" << std::endl;
     exit(0);
   }

};


void GreensTensorVacuum::integrate_k_1d(cx_mat::fixed<3,3>& GT, double omega, Options opts)
{
  GT.zeros();
  opts.omega = omega;
  ////opts.class_pt = this;

  if (opts.fancy_R)
  {
    std::cerr<< "Real part of the vacuum Green's is divergent and already included within the omega_a" << std::endl;
    exit(0);
  }

  opts.indices(0) = 0;
  opts.indices(1) = 0;
  GT(0,0) = cquad(&integrand_k_1d,&opts,-omega/(1.0+this->v),omega/(1.0-this->v),1E-7,0);
  opts.indices(0) = 1;
  opts.indices(1) = 1;
  GT(1,1) = cquad(&integrand_k_1d,&opts, -omega/(1.0+this->v),omega/(1.0-this->v),1E-7,0);
  GT(2,2) = GT(1,1);
};

double GreensTensorVacuum::integrand_k_1d(double kv, void *opts)
{
  //Units: c=1, 4 pi epsilon_0 = 1, hbar = 1
  Options* opts_pt = static_cast<Options*>(opts);
  GreensTensorVacuum* pt = static_cast<GreensTensorVacuum*>(opts_pt->class_pt);
  double omega = opts_pt->omega;
  double result;
  double xi_quad, omega_pl, omega_pl_quad;
  double beta = pt->beta;
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
