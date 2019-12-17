// ini parser
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
namespace pt = boost::property_tree;

#include "PolarizabilityBath.h"
#include "MemoryKernel/MemoryKernelFactory.h"

PolarizabilityBath::PolarizabilityBath(double omega_a, double alpha_zero, MemoryKernel *mu, GreensTensor *green)
{
  // set parameters
  this->omega_a = omega_a;
  this->alpha_zero = alpha_zero;

  // set memory kernel
  this->mu = mu;

  // set green's function
  this->greens_tensor = green;
}



PolarizabilityBath::PolarizabilityBath(std::string input_file)
{
  // Create a root
  pt::ptree root;

  // Load the ini file in this ptree
  pt::read_ini(input_file, root);

  // read parameters
  this->omega_a = root.get<double>("Polarizability.omega_a");
  this->alpha_zero = root.get<double>("Polarizability.alpha_zero");

  // read memory kernel
  this->mu = MemoryKernelFactory::create(input_file);
};


std::complex<double> PolarizabilityBath::get_mu(double omega)
{
    return this->mu->mu(omega);
};



void PolarizabilityBath::calculate(cx_mat::fixed<3,3>& alpha, double omega)
{
  // imaginary unit
  std::complex<double> I(0.0, 1.0);

  // calculate diagonal entries
  cx_mat::fixed<3,3> diag;
  diag.zeros();
  diag(0,0) = diag(1,1) = diag(2,2) = omega_a*omega_a - omega*omega - I * omega * mu->mu(omega);

  // calculate integral over green's tensor with fancy R
  cx_mat::fixed<3,3> greens_R;
  struct Options_GreensTensor opts_R;
  opts_R.fancy_R = true;
  opts_R.omega = omega;
  opts_R.class_pt = this->greens_tensor;
  this->greens_tensor->integrate_k_1d(greens_R, opts_R);

  // calculate integral over green's tensor with fancy I
  cx_mat::fixed<3,3> greens_I;
  struct Options_GreensTensor opts_I;
  opts_I.fancy_I = true;
  opts_I.omega = omega;
  opts_I.class_pt = this->greens_tensor;

  this->greens_tensor->integrate_k_1d(greens_I, opts_I);

  // put everything together
  alpha = alpha_zero * omega_a * omega_a * inv(diag - alpha_zero * omega_a * omega_a * (greens_R + I*greens_I));
};
