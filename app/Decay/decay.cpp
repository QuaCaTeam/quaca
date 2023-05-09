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
        "threads", po::value<int>(&(num_threads))->default_value(1),
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

  // read parameters (No Looper class used for omega, since it is no attribute of the class)
  
  pt::ptree root;
  pt::read_json(parameter_file, root);
  
  double start;              // starting value
  double end;                // end value
  int number_of_steps;       // number of steps
  std::string scale;         // scale type

  start = root.get<double>("Looper.start");
  end = root.get<double>("Looper.end");
  number_of_steps = root.get<double>("Looper.steps");
  scale = root.get<std::string>("Looper.scale");
  
  // Calculate steps
  std::vector<double> steps; // array containing the steps
  if (scale == "linear") {
    double spacing = (end - start) / ((double)number_of_steps - 1);
    for (int i = 0; i < number_of_steps; i++) {
      steps.push_back(start + i * spacing);
    }
  } else if (scale == "log") {
    double spacing = pow(end / start, 1. / ((double)number_of_steps - 1.0));
    for (int i = 0; i < number_of_steps; i++) {
      steps.push_back(start * pow(spacing, i));
    }

  } else {
    std::cerr << "Unknown scale: " << scale << std::endl;
    exit(-1);
  }


  // array to store all the computed values of the loop
  std::vector<double> decay_data;
  decay_data.resize(number_of_steps);

  // define progressbar
  ProgressBar progbar(number_of_steps, 70);

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

    // define needed quantities
    auto polarizability = std::make_shared<Polarizability>(parameter_file);
    double omega_a = polarizability->get_omega_a();
    double alpha_zero = polarizability->get_alpha_zero();
    // Parallelize the for-loop of the given looper
#pragma omp critical
    progbar.display();

#pragma omp for schedule(dynamic)
    for (int i = 0; i < number_of_steps; i++) {

	    // Calculate decay rate   
	    // define different alphas
	    std::complex<double> I(0e0, 1e0);
	    cx_mat::fixed<3, 3> alphaI;
	    cx_mat::fixed<3, 3> alphaR;
	    cx_mat::fixed<3, 3> inv_alpha;
	    cx_mat::fixed<3, 3> inv_alpha_dag;

	    // calculate alpha
	    polarizability->calculate_tensor(steps[i], alphaI, IM);
	    polarizability->calculate_tensor(steps[i], alphaR, RE);

	    inv_alpha     = inv(alphaR + I * alphaI);
	    inv_alpha_dag = inv(alphaR - I * alphaI);

	    decay_data[i] = alpha_zero*omega_a*omega_a*real(trace(inv_alpha*alphaI*inv_alpha_dag))/steps[i];

#pragma omp critical
	    ++progbar;
#pragma omp critical
	    progbar.display();
    }
  }

  // define output file
  std::ofstream file;
  file.open(output_file);

  // write the data to the file
  double step;
  for (int i = 0; i < number_of_steps; i++) {
    step = steps[i];
    file << step << "," << decay_data[i] << "\n";
  }

  // close file
  file.close();

  // close progress bar
  progbar.done();

  return 0;
}
