// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "Looper.h"
#include <cassert>

Looper::Looper(double start, double end, int number_of_steps, const std::string &scale)
    : start(start), end(end), number_of_steps(number_of_steps), scale(scale) {
  assert(start < end);
  this->calculate_steps();
}

Looper::Looper(const std::string &input_file) {

  // read parameters
  pt::ptree root;
  pt::read_json(input_file, root);
  this->start = root.get<double>("Looper.start");
  this->end = root.get<double>("Looper.end");
  this->number_of_steps = root.get<double>("Looper.steps");
  this->scale = root.get<std::string>("Looper.scale");

  assert(start < end);
  this->calculate_steps();
}

void Looper::calculate_steps() {
  if (scale == "linear") {
    double spacing = (end - start) / ((double)number_of_steps - 1);
    for (int i = 0; i < number_of_steps; i++) {
      this->steps.push_back(start + i * spacing);
    }
  } else if (scale == "log") {
    double spacing = pow(end / start, 1. / ((double)number_of_steps - 1.0));
    for (int i = 0; i < number_of_steps; i++) {
      this->steps.push_back(start * pow(spacing, i));
    }

  } else {
    std::cerr << "Unknown scale: " << scale << std::endl;
    exit(-1);
  }
}
