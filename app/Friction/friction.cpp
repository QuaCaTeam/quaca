#include <iostream>
#include <omp.h>

// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "ProgressBar.hpp"
#include "Quaca.h"

int main(int argc, char *argv[]) {
  // get file from command line
  Options opts(argc, argv);
  std::string parameters = opts.get_parameter_file();

  // define looper
  auto looper = LooperFactory::create(parameters);

  // array to store all the computed values of the loop
  std::vector<double> friction_data;
  friction_data.resize(looper->get_steps_total());

  // define progressbar
  ProgressBar progbar(looper->get_steps_total(), 70);

  // Check wether the number of threads have been set by the --threads flag
  if (opts.get_num_threads() == -1) {
    // If the --threads flag has not been set, set the number of flags to
    // the maximal value
    opts.set_num_threads(omp_get_max_threads());
  }
  // Create a paralle region given threads given by the --threads flag
  // we have to create the parallel region already here to ensure,
  // that any thread creates their own instance of quantum_friction
  std::cout << "Starting parallel region with " << opts.get_num_threads()
            << " threads." << std::endl;
  if (opts.get_num_threads() > omp_get_max_threads()) {
    std::cout
        << "There are not enough avaiable threads. Maximal avaiable threads: "
        << omp_get_max_threads() << std::endl;
    std::cout << "Aborting calculation" << std::endl;
    exit(0);
  }
#pragma omp parallel num_threads(opts.get_num_threads())
  {

    // Create a root
    pt::ptree root;

    // Load the json file in this ptree
    pt::read_json(parameters, root);

    double relerr_omega = root.get<double>("Friction.relerr_omega");

    // define needed quantities
    auto polarizability = std::make_shared<Polarizability>(parameters);
    auto powerspectrum = std::make_shared<PowerSpectrum>(
        polarizability->get_greens_tensor(), polarizability);
    auto quant_friction =
        std::make_shared<Friction>(polarizability->get_greens_tensor(),
                                   polarizability, powerspectrum, relerr_omega);

    // Parallelize the for-loop of the given looper
    progbar.display();

#pragma omp for schedule(dynamic)
    for (int i = 0; i < looper->get_steps_total(); i++) {
      friction_data[i] = looper->calculate_value(i, quant_friction);

      // std::cout << "Thread: " << omp_get_thread_num() << " Step: " << i
      //           << " Friction force: " << friction_data[i] << endl;
      ++progbar;
      progbar.display();
    }
  }

  // Store calculated data

  // define output file
  std::ofstream file;
  file.open(opts.get_output_file());

  // write the data to the file
  double step;
  for (int i = 0; i < looper->get_steps_total(); i++) {
    step = looper->get_step(i);
    file << step << "," << friction_data[i] << "\n";
  }

  // close file
  file.close();

  // close progress bar
  progbar.done();

  return 0;
}
