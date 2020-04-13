#include "PowerSpectrum.h"
#include "../GreensTensor/GreensTensorFactory.h"
#include "../Polarizability/PolarizabilityFactory.h"
#include <string>

// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <utility>
namespace pt = boost::property_tree;

// Constructor using a json-file
PowerSpectrum::PowerSpectrum(const std::string &input_file) {
  // Create a root
  pt::ptree root;

  // Load the json file in this ptree
  pt::read_json(input_file, root);

  // initialize Green's tensor by an input file
  this->greens_tensor = GreensTensorFactory::create(input_file);
  // initialize polarizability by an input file
  this->polarizability = PolarizabilityFactory::create(input_file);
}

// Constructor with initialization list
PowerSpectrum::PowerSpectrum(std::shared_ptr<GreensTensor> greens_tensor,
                             std::shared_ptr<Polarizability> polarizability)
    : greens_tensor(std::move(greens_tensor)),
      polarizability(std::move(polarizability)){}
