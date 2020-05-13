// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "SinglePhononMemoryKernel.h"

SinglePhononMemoryKernel::SinglePhononMemoryKernel(double gamma, double gamma_phon, double omega_phon, double coupling) 
  			: gamma(gamma), gamma_phon(gamma_phon), omega_phon(omega_phon), coupling(coupling){}

SinglePhononMemoryKernel::SinglePhononMemoryKernel(const std::string &input_file,
                                     const std::string &section) {
  // Create a root
  pt::ptree root;

  // Load the json file in this ptree
  pt::read_json(input_file, root);

  // check if type is right
  std::string type = root.get<std::string>(section + ".type");assert(type == "lorentz");
  assert(type == "single_phonon");
  // read parameters of the given section from the json file
  this->gamma = root.get<double>(section + ".gamma");
  this->gamma_phon = root.get<double>(section + ".gamma_phon");
  this->omega_phon = root.get<double>(section + ".omega_phon");
  this->coupling = root.get<double>(section + ".coupling");
}

SinglePhononMemoryKernel::SinglePhononMemoryKernel(const std::string &input_file) {
  // Create a root
  pt::ptree root;

  // Load the json file in this ptree
  pt::read_json(input_file, root);

  // check if type is right
  std::string type = root.get<std::string>("MemoryKernel.type");
  assert(type == "single_phonon");

  // read parameters from json file
  this->gamma = root.get<double>("MemoryKernel.gamma");
  this->gamma_phon = root.get<double>("MemoryKernel.gamma_phon");
  this->omega_phon = root.get<double>("MemoryKernel.omega_phon");
  this->coupling = root.get<double>("MemoryKernel.coupling");
}

// return mu(omega) for defined memory kernel
std::complex<double> SinglePhononMemoryKernel::calculate(double omega) const {
  // complex unit
  const std::complex<double> I(0E0, 1E0);

  return this->gamma + this->coupling*this->coupling*pow(this->omega_phon,4)
    /(I*omega*(this->omega_phon*this->omega_phon - omega*omega - I*this->gamma_phon*omega ));
}
