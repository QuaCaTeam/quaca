#include "GreensTensorVacuum.h"

GreensTensorVacuum::GreensTensorVacuum(double a)
{
  // set velocity
  this->v = a;

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
void GreensTensorVacuum::calculate_integrated(cx_mat::fixed<3,3>& GT, double omega, double kv, std::vector<std::string> options)
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


void GreensTensorVacuum::calculate_integrated(cx_mat::fixed<3,3>& GT, double omega, std::vector<std::string> options)
{
 // ONLY FOR TEST PURPOSE
  GreensTensorVacuum::calculate_pure(GT, vec(2,fill::zeros) ,omega);
};
