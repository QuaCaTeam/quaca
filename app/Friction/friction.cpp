#include <iostream>
#include <omp.h>

// program options
#include <boost/program_options.hpp>
namespace po = boost::program_options;

// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "ProgressBar.hpp"
#include "Quaca.h"

// parameters that are parsed from the command line
std::string parameter_file;
std::string output_file;
int num_threads = 1;

// reads the input file and number of threads from the command line
// uses boost program options
void read_command_line(int argc, char *argv[]) {
  /* Read command line options */
  try {
    // List all options and their description
    po::options_description desc("Allowed options");
    desc.add_options()("help,h", "Help screen")(
        "file", po::value<std::string>(&(parameter_file)), "Input File")(
        "threads", po::value<int>(&(num_threads))->default_value(-1),
        "Number of parallel threads");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // if the help option is given, show the flag description
    if (vm.count("help")) {
      std::cout << desc << std::endl;
      exit(0);
    }

  } catch (std::exception &e) {
    std::cerr << "error: " << e.what() << std::endl;
  }

  /* set according output file */
  output_file =
      parameter_file.substr(0, parameter_file.find_last_of('.')) +
      ".csv";
}

int main(int argc, char *argv[]) {
  // get command line options
  read_command_line(argc, argv);

  // define looper
  auto looper = LooperFactory::create(parameter_file);

  // array to store all the computed values of the loop
  std::vector<double> friction_data;
  friction_data.resize(looper->get_steps_total());

  // define progressbar
  ProgressBar progbar(looper->get_steps_total(), 70);


  // Check whether the number of threads have been set by the --threads flag
  if (num_threads == -1) {
    // If the --threads flag has not been set, set the number of flags to
    // the maximal value
    num_threads = omp_get_max_threads();
  }

  // Create a parallel region given threads given by the --threads flag
  // we have to create the parallel region already here to ensure,
  // that any thread creates their own instance of quantum_friction
  std::cout << "Starting parallel region with " << num_threads
            << " threads." << std::endl;
  if (num_threads > omp_get_max_threads()) {
    std::cout
        << "There are not enough avaiable threads. Maximal avaiable threads: "
        << omp_get_max_threads() << std::endl;
    std::cout << "Aborting calculation" << std::endl;
    exit(0);
  }
#pragma omp parallel num_threads(num_threads)
  {

    // Create a root
    pt::ptree root;

    // Load the json file in this ptree
    pt::read_json(parameter_file, root);

    double relerr_omega = root.get<double>("Friction.relerr_omega");

    // define needed quantities
    auto polarizability = std::make_shared<Polarizability>(parameter_file);
    auto powerspectrum = std::make_shared<PowerSpectrum>(
        polarizability->get_greens_tensor(), polarizability);
    auto quant_friction =
        std::make_shared<Friction>(polarizability->get_greens_tensor(),
                                   polarizability, powerspectrum, relerr_omega);

    // Parallelize the for-loop of the given looper
#pragma omp critical
    progbar.display();

#pragma omp for schedule(dynamic)
    for (int i = 0; i < looper->get_steps_total(); i++) {
      friction_data[i] = looper->calculate_value(i, quant_friction);
#pragma omp critical
      {
          // define output file
          std::ofstream file;
          file.open(output_file);

          // write results in output file
          for(int j = 0; j < looper->get_steps_total(); j++) {
            double step = looper->get_step(j);
            file << step << "," << friction_data[j] << "\n";
          }
          // close file
          file.close();
      }

#pragma omp critical
      ++progbar;
#pragma omp critical
      progbar.display();
    }
  }


  // close progress bar
  progbar.done();

  return 0;
}
