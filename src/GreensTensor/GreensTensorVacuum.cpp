#include "GreensTensorVacuum.h"
#include "../Calculations/Integrations.h"

GreensTensorVacuum::GreensTensorVacuum(double v, double beta)
{
  // set velocity
  this->v = v;
  // set inverse temperature
  this->beta = beta;

};
void GreensTensorVacuum::calculate_pure(cx_mat::fixed<3,3>& GT, vec::fixed<2> kvec, double omega)
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
void GreensTensorVacuum::calculate_integrated(cx_mat::fixed<3,3>& GT, double omega, double kv)
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
};


void GreensTensorVacuum::calculate_integrated(cx_mat::fixed<3,3>& GT, double omega, Options *opts)
{
  // If false the G_{xx} component is calculated.
  // For true the G_{yy/zz} component is calculated.
  bool zz_component;

  GT.zeros();
  
  if (opts->fancy_R)
  {
    std::cerr<< "Real part of the vacuum Green's is divergent and already included within the omega_a" << std::endl;
    exit(0);
  };

  double fint(double kv, void *p)
  {
    double result;
    double xi2, omegapl, omegapl2;
    omegapl = (omega+kv*this->v);
    omegapl2 = omegapl*omegapl;
    xi2 = omegapl2 - kv*kv;

    if (zz_component)
    {
      result = 0.5*(omegapl2 - xi2*0.5);
    };
    else
    {
      result = 0.5*xi2;
    };
    if (fancy_I)
    {
      return result;
    };
    if (fancy_I_kv)
    {
      result = result*kv;
    }
    else if (fancy_I_temp)
    {
      result = result/(1.0-exp(-beta*omegapl))
    }
    else if (fancy_I_kv_temp)
    {
      result = result*kv/(1.0-exp(-beta*omegapl))
    };
    return result;
  };

  zz_component = false;
  GT(0,0) = cquad(fint,-omega/(1.0+v),omega/(1.0-v),1E-7,0);
  zz_component = true;
  GT(1,1) = cquad(fint,-omega/(1.0+v),omega/(1.0-v),1E-7,0);
  GT(2,2) = GT(1,1);
};
