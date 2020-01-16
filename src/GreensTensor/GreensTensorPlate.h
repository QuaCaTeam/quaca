#ifndef GREENSTENSORPLATE_H
#define GREENSTENSORPLATE_H

#include <complex>
#include <cmath>
#include "GreensTensor.h"
#include "Permittivity/PermittivityFactory.h"
// ini parser
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
namespace pt = boost::property_tree;

class GreensTensorPlate : public GreensTensor
{
private:
    // permittivity is needed to describe the surface's response
    Permittivity *permittivity;
    // kappa_cut defines the numerical cut-off of the kappa integration
    double delta_cut;
    double za;
public:
  GreensTensorPlate(double v, double za, double beta, Permittivity *permittivity, double delta_cut):
  GreensTensor(v, beta) {this->za = za; this->delta_cut = delta_cut; this->permittivity = permittivity;};


  GreensTensorPlate(std::string input_file):
  GreensTensor(input_file){
  this->permittivity = PermittivityFactory::create(input_file);

  // Create a root
  pt::ptree root;

  // Load the ini file in this ptree
  pt::read_ini(input_file, root);

  // read parameters
  this->za = root.get<double>("GreensTensor.za");
  this->delta_cut = root.get<double>("GreensTensor.delta_cut");
  };

  void calculate_tensor(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts);
  void integrate_k_2d(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts);
  void integrate_k_1d(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts);

  std::complex<double> get_epsilon(double omega){return this->permittivity->epsilon(omega);};
  double get_za(){return this->za;}
  double get_delta_cut(){return this->delta_cut;}
  static double integrand_k_1d(double kx, void* opts);
  static double integrand_k_2d(double ky, void* opts);

};


#endif // GREENSTENSORPLATE_H
