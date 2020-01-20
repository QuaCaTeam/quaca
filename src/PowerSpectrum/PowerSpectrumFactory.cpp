#include <iostream>

// ini parser
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "PowerSpectrumFactory.h"
#include "PowerSpectrumHarmOsc.h"

// power spectrum
PowerSpectrum *PowerSpectrumFactory::create(std::string input_file) {
  // set return pointer to NULL
  PowerSpectrum *powerspectrum = NULL;

  // Create a root
  pt::ptree root;

  // Load the ini file in this ptree
  pt::read_ini(input_file, root);

  // read the type of the power spectrum
  std::string type = root.get<std::string>("PowerSpectrum.type");

  // set the right pointer, show error if type is unknown
  if (type == "harmonic oscillator") {
    powerspectrum = new PowerSpectrumHarmOsc(input_file);
  }

  else {
    std::cerr << "Error: Unknown power spectrum type (" << type << ")!"
              << std::endl;
    exit(0);
  };

  // return powerspectrum pointer
  return powerspectrum;
};
