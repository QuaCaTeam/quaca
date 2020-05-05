// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "LorentzMemoryKernel.h"

LorentzMemoryKernel::LorentzMemoryKernel(double gamma, double omega_0, double omega_p, double eps_inf) 
  			: gamma(gamma), omega_0(omega_0), omega_p(omega_p), eps_inf(eps_inf){}

LorentzMemoryKernel::LorentzMemoryKernel(const std::string &input_file,
                                     const std::string &section) {
  // Create a root
  pt::ptree root;

  // Load the json file in this ptree
  pt::read_json(input_file, root);

  // check if type is right
  std::string type = root.get<std::string>(section + ".type");
  assert(type == "lorentz");

  // read parameters of the given section from the json file
  this->gamma = root.get<double>(section + ".gamma");
  this->omega_0 = root.get<double>(section + ".omega_0");
  this->omega_p = root.get<double>(section + ".omega_p");
  this->eps_inf = root.get<double>(section + ".eps_inf");
}

LorentzMemoryKernel::LorentzMemoryKernel(const std::string &input_file) {
  // Create a root
  pt::ptree root;

  // Load the json file in this ptree
  pt::read_json(input_file, root);

  // check if type is right
  std::string type = root.get<std::string>("MemoryKernel.type");
  assert(type == "lorentz");

  // read parameters from json file
  this->gamma = root.get<double>("MemoryKernel.gamma");
  this->omega_0 = root.get<double>("MemoryKernel.omega_0");
  this->omega_p = root.get<double>("MemoryKernel.omega_p");
  this->eps_inf = root.get<double>("MemoryKernel.eps_inf");
}

// return mu(omega) for defined memory kernel
std::complex<double> LorentzMemoryKernel::calculate(double omega) const {
  // complex unit
  const std::complex<double> I(0E0, 1E0);

  return this->eps_inf + this->omega_p*this->omega_p
    /(this->omega_0*this->omega_0 - omega*omega +I*this->gamma*omega );
}
