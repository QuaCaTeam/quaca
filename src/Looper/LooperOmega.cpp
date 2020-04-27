// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "LooperOmega.h"
#include "../Polarizability/Polarizability.h"

LooperOmega::LooperOmega(double start, double end, int number_of_steps,
                 std::string scale)
    : Looper(start, end, number_of_steps, scale){};

LooperOmega::LooperOmega(std::string input_file) : Looper(input_file) {
  // Create a root
  pt::ptree root;

  // Load the json file in this ptree
  pt::read_json(input_file, root);

  // check if type is right
  std::string type = root.get<std::string>("Looper.type");
  assert(type == "omega");
};

double LooperOmega::calculate_value(int step, void *quantity) {

  Polarizability *polarizability = static_cast<Polarizability *>(quantity);

    // imaginary unit
  std::complex<double> I(0.0, 1.0);

  // define different alphas
  cx_mat::fixed<3, 3> alphaI;
  cx_mat::fixed<3, 3> alphaR;
  cx_mat::fixed<3, 3> inv_alpha;
  cx_mat::fixed<3, 3> inv_alpha_dag;

  // set the option struct
  Options_Polarizability opts;
  opts.class_pt = polarizability;

  // change omega
  opts.omega = this->steps[step];

  // calculate alpha
  opts.fancy_complex = IM;
  polarizability->calculate_tensor(alphaI, opts);
  opts.fancy_complex = RE;
  polarizability->calculate_tensor(alphaR, opts);

  inv_alpha     = inv(alphaR + I * alphaI);
  inv_alpha_dag = inv(alphaR - I * alphaI);

  // calculate the trace


  return real(trace(inv_alpha*alphaI*inv_alpha_dag))/(opts.omega*3);
};

