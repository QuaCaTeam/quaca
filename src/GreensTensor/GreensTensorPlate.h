#ifndef GREENSTENSORPLATE_H
#define GREENSTENSORPLATE_H

// ini parser
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
namespace pt = boost::property_tree;

#include "GreensTensor.h"

class GreensTensorPlate : public GreensTensor
{
private:
  double z_a; // distance from plate

public:

  // constructors
  GreensTensorPlate(std::string input_file): GreensTensor(input_file)
  {
    pt::ptree root;
    pt::read_ini(input_file, root);
    this->z_a = root.get<double>("GreensTensor.z_a");
  };

  GreensTensorPlate(double v, double z_a, double beta): GreensTensor(v, beta), z_a(z_a)
  {
    assert(z_a > 0);
  };

  // calculate full tensor
  void calculate_tensor(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts);

  // integrate over a two-dimensional k space
  void integrate_2d_k(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts);

  // integrate over a one-dimensional k space
  void integrate_1d_k(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts);

};


#endif // GREENSTENSORPLATE_H
