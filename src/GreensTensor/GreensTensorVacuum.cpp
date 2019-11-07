#include "GreensTensorVacuum.h"

GreensTensorVacuum::GreensTensorVacuum(double a)
{
  // set velocity
  this->v = a;

};
void GreensTensorVacuum::calculate_pure(cx_mat::fixed<3,3>& GT, vec::fixed<2> kvec, double omega)
{
  // calculating the solely the imaginary part of the free Green tensor
  std::complex<double> diag;
  diag = omega*omega*omega*2.0/3.0;
  GT.zeros();
  GT(0,0) = diag;
  GT(1,1) = diag;
  GT(2,2) = diag;
};
void GreensTensorVacuum::calculate_integrated(cx_mat::fixed<3,3>& GT, double omega, std::vector<std::string> options)
{
 // ONLY FOR TEST PURPOSE
  GreensTensorVacuum::calculate_pure(GT, vec(2,fill::zeros) ,omega);
};
