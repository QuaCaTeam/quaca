#include <cassert>

// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "../Calculations/Integrations.h"
#include "GreensTensorVacuum.h"

GreensTensorVacuum::GreensTensorVacuum(double v, double beta, double relerr)
    : GreensTensor(v, beta), relerr(relerr) {}

GreensTensorVacuum::GreensTensorVacuum(const std::string &input_file)
    : GreensTensor(input_file) {

  // Create a root
  pt::ptree root;

  // Load the json file in this ptree
  pt::read_json(input_file, root);

  // Load relative accuracy
  this->relerr = root.get<double>("GreensTensor.rel_err_1");

  // check if type is right
  std::string type = root.get<std::string>("GreensTensor.type");
  assert(type == "vacuum");
}

// Compute the full Green's tensor for a given frequency \omega and a given
// momentum vector k For the definition see notes/VacuumGreen.pdf eq. (2)
void GreensTensorVacuum::calculate_tensor(cx_mat::fixed<3, 3> &GT,
                                          Options_GreensTensor opts) {
  if (opts.fancy_complex == IM) {
    // Read out the k-vector and the frequency omega
    double k_x = opts.kvec(0);
    double k_y = opts.kvec(1);
    double omega = opts.omega;

    // Define useful variables
    double k_quad = k_x * k_x + k_y * k_y;
    double omega_quad = omega * omega;

    // Reset tensor in which the final result is stored
    GT.zeros();

    // Ensure that heavyside function is fulfilled
    if (omega_quad - k_quad > 0) {
      // Compute the diagonal components of the tensor. The off-diagonal
      // elements are all zero
      double pre = 1.0 / (2 * M_PI * sqrt(omega_quad - k_quad));
      GT(0, 0) = pre * (omega_quad - k_x * k_x);
      GT(1, 1) = pre * (omega_quad - k_y * k_y);
      GT(2, 2) = pre * k_quad;
    }
  } else {
    std::cerr << "Only the imaginary part of the Green tensor with Doppler "
                 "shift in the frequency argument is implemented"
              << std::endl;
    exit(0);
  }
}

// Compute the integration with respect to the 2-d k vector
// Ref: notes/VacuumFriction.pdf eq. (10)
void GreensTensorVacuum::integrate_k(cx_mat::fixed<3, 3> &GT,
                                     Options_GreensTensor opts) {
  if (opts.fancy_complex == RE) {
    // Even though the real part of the Green's tensor is not implemented, a
    // default return value of an empty tensor was chosen, to allow for the
    // general structure of the polarizability to depend both on the real and
    // imaginary part of a given Green's tensor
    GT.zeros();
  } else if (opts.fancy_complex == IM) {
    // Read out the frequency \omega
    double omega = opts.omega;

    // Reset the tensor to store the final result
    GT.zeros();
    // Ensure that the integration limits are properly ordered
    if (omega >= 0) {

      // Numerically integrate the xx component
      opts.indices(0) = 0;
      opts.indices(1) = 0;
      GT(0, 0) = cquad(&integrand_k, &opts, -omega / (1.0 + this->v),
                       omega / (1.0 - this->v), this->relerr, 0);
      opts.indices(0) = 1;
      opts.indices(1) = 1;

      GT(1, 1) = cquad(&integrand_k, &opts, -omega / (1.0 + this->v),
                       omega / (1.0 - this->v), this->relerr, 0);
      GT(2, 2) = GT(1, 1);
    }
    // Switching the integration bounds for negative frequencies
    if (omega < 0) {

      // Numerically integrate the xx component
      opts.indices(0) = 0;
      opts.indices(1) = 0;
      GT(0, 0) = -cquad(&integrand_k, &opts, omega / (1.0 - this->v),
                        -omega / (1.0 + this->v), this->relerr, 0);
      opts.indices(0) = 1;
      opts.indices(1) = 1;

      GT(1, 1) = -cquad(&integrand_k, &opts, omega / (1.0 - this->v),
                        -omega / (1.0 + this->v), this->relerr, 0);
      GT(2, 2) = GT(1, 1);
    }
  }
}

// Implementation of the different integrands for the integration
// of the 2-d k-vector
// Ref: notes/VacuumFriction eq. (10) and (11)
double GreensTensorVacuum::integrand_k(double kv, void *opts) {
  // Casting the class-pointer to the correct pointer-type
  auto *opts_pt = static_cast<Options_GreensTensor *>(opts);
  auto *pt = dynamic_cast<GreensTensorVacuum *>(opts_pt->class_pt);

  // Read out the relevant parameters
  double omega = opts_pt->omega;
  double beta = pt->beta;

  double omega_pl = (omega + kv * pt->v);
  double omega_pl_quad = omega_pl * omega_pl;
  double xi_quad = omega_pl_quad - kv * kv;

  // Variable to store the final result
  double result = 0;

  // Only the imaginary part is implemented
  if (opts_pt->fancy_complex == IM) {
    // Compute the basis integrand of eq. (10)
    if (opts_pt->indices(0) == 0 && opts_pt->indices(1) == 0) {
      result = 0.5 * xi_quad;
    } else if (opts_pt->indices(0) == opts_pt->indices(1)) {
      result = 0.5 * (omega_pl_quad - xi_quad * 0.5);
    } else {
      return 0;
    }

    // Multply with the additional weight function f, the options can be found
    // in eq. (11)
    if (opts_pt->weight_function == KV) {
      result *= kv;
    } else if (opts_pt->weight_function == TEMP) {
      result /= (1.0 - exp(-beta * omega_pl));
    } else if (opts_pt->weight_function == KV_TEMP) {
      result *= kv / (1.0 - exp(-beta * omega_pl));
    } else if (opts_pt->weight_function == NON_LTE) {
      result *=
          (1. / (1. - exp(-beta * omega_pl)) - 1. / (1. - exp(-beta * omega)));
    } else if (opts_pt->weight_function == KV_NON_LTE) {
      result *= kv * (1. / (1. - exp(-beta * omega_pl)) -
                      1. / (1. - exp(-beta * omega)));
    }
  }

  return result;
}

double GreensTensorVacuum::omega_ch() {
  double result = 0;
  return result;
}

void GreensTensorVacuum::print_info(std::ofstream &file) {
  file << "# GreensTensorVacuum "
       << "\n"
       << "# v = " << v << "\n"
       << "# beta = " << beta << "\n"
       << "# rel_err_1 = " << relerr << "\n"
       << "\n";
}
