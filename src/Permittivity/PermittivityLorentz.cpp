// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <utility>
namespace pt = boost::property_tree;

#include "../MemoryKernel/MemoryKernelFactory.h"
#include "PermittivityLorentz.h"

PermittivityLorentz::PermittivityLorentz(
    double eps_inf, double omega_p, double omega_0,
    std::shared_ptr<MemoryKernel> memory_kernel)
    : eps_inf(eps_inf), omega_p(omega_p), omega_0(omega_0),
      memory_kernel(std::move(memory_kernel)) {}

// constructor for drude model from .json file
PermittivityLorentz::PermittivityLorentz(const std::string &input_file) {

  // Create a root
  pt::ptree root;

  // Load the json file in this ptree
  pt::read_json(input_file, root);

  // check if type is right
  std::string type = root.get<std::string>("Permittivity.type");
  assert(type == "lorentz");

  // read parameters
  this->eps_inf = root.get<double>("Permittivity.eps_inf");
  this->omega_p = root.get<double>("Permittivity.omega_p");
  this->omega_0 = root.get<double>("Permittivity.omega_0");

  this->memory_kernel =
      MemoryKernelFactory::create(input_file, "Permittivity.MemoryKernel");
}

// calculate the permittivity
std::complex<double> PermittivityLorentz::calculate(double omega) const {
  // dummies for result and complex unit
  std::complex<double> result;
  std::complex<double> I(0.0, 1.0);

  // calculate the result
  result = eps_inf - omega_p * omega_p /
                         (omega_0 * omega_0 - omega * omega -
                          I * omega * memory_kernel->calculate(omega));

  return result;
}

// calculate the permittivity scaled by omega
std::complex<double>
PermittivityLorentz::calculate_times_omega(double omega) const {
  // dummies for result and complex unit
  std::complex<double> result;
  std::complex<double> I(0.0, 1.0);

  // calculate the result
  result = eps_inf * omega - omega_p * omega_p * omega /
                                 (omega_0 * omega_0 - omega * omega -
                                  I * omega * memory_kernel->calculate(omega));

  return result;
}

void PermittivityLorentz::print_info(std::ostream &stream) const {
  stream << "# PermittivityLorentz\n#\n"
         << "# eps_inf = " << eps_inf << "\n"
         << "# omega_p = " << omega_p << "\n"
         << "# omega_0 = " << omega_0 << "\n";
  memory_kernel->print_info(stream);
}
