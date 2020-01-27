// ini parser
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "../GreensTensor/GreensTensorFactory.h"
#include "Polarizability.h"

Polarizability::Polarizability(double omega_a, double alpha_zero,
                               GreensTensor *greens_tensor)
    : omega_a(omega_a), alpha_zero(alpha_zero), greens_tensor(greens_tensor){};

Polarizability::Polarizability(std::string input_file) {
  // Create a root
  pt::ptree root;

  // Load the ini file in this ptree
  pt::read_ini(input_file, root);

  // read parameters
  this->omega_a = root.get<double>("Polarizability.omega_a");
  this->alpha_zero = root.get<double>("Polarizability.alpha_zero");

  // read greens tensor
  this->greens_tensor = GreensTensorFactory::create(input_file);

  // read type
  this->type = root.get<std::string>("Polarizability.type");
};

double Polarizability::integrand_omega(double omega, void *opts) {
  Options_Polarizability *opts_pt = static_cast<Options_Polarizability *>(opts);
  cx_mat::fixed<3, 3> alpha;

  opts_pt->omega = omega;
  opts_pt->class_pt->calculate_tensor(alpha, *opts_pt);

  assert(imag(alpha(opts_pt->indices(0), opts_pt->indices(1))) == 0);
  double result = real(alpha(opts_pt->indices(0), opts_pt->indices(1)));

  return result;
};

double Polarizability::integrate_omega(Options_Polarizability opts,
                                       double omega_min, double omega_max,
                                       double relerr, double abserr) {
  return cquad(&integrand_omega, &opts, omega_min, omega_max, relerr, abserr);
}
