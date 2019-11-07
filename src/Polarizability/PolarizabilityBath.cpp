// ini parser
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
namespace pt = boost::property_tree;

#include "PolarizabilityBath.h"
#include "MemoryKernel/MemoryKernelFactory.h"

PolarizabilityBath::PolarizabilityBath(double a, double b, MemoryKernel *mu)
{
  // set parameters
  this->omega_a = a;
  this->alpha_zero = b;

  // set memory kernel
  this->memorykernel = mu;
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
  this->memorykernel = MemoryKernelFactory::create(input_file);
};

std::complex<double> PolarizabilityBath::get_mu(double omega)
{
    return this->memorykernel->mu(omega);
};


void PolarizabilityBath::calculate(cx_mat::fixed<3,3>& alpha, double omega)
{
  // calculate entries
  std::complex<double> diag;
  std::complex<double> I(0.0, 1.0);
  diag = 1.0/ (omega_a*omega_a - omega*omega - I * omega * memorykernel->mu(omega));

  // set alpha matrix entries
  alpha(0,0) = diag;
  alpha(1,1) = diag;
  alpha(2,2) = diag;
};
