#include <assert.h>

// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "../Calculations/Integrations.h"
#include "GreensTensorVacuum.h"

GreensTensorVacuum::GreensTensorVacuum(double v, double beta, double relerr)
    : GreensTensor(v, beta), relerr(relerr) {
    assert(relerr >= 0);
    }

GreensTensorVacuum::GreensTensorVacuum(std::string input_file)
    : GreensTensor(input_file) {

  // Create a root
  pt::ptree root;

  // Load the json file in this ptree
  pt::read_json(input_file, root);

  // Load relative accuracy
  this->relerr = root.get<double>("GreensTensor.rel_err_1");

  assert(relerr >= 0);

  // check if type is right
  std::string type = root.get<std::string>("GreensTensor.type");
  assert(type == "vacuum");
}

// Compute the full Green's tensor for a given frequency \omega and a given
// momentum vector k For the definition see notes/VacuumGreen.pdf eq. (2)
void GreensTensorVacuum::calculate_tensor(double omega, vec::fixed<2> k,
                                          cx_mat::fixed<3, 3> &GT) const {
  // Read out the k-vector and the frequency \omega
  double k_x = k(0);
  double k_y = k(1);

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
}

// Compute the integration with respect to the 2-d k vector
// Ref: notes/VacuumFriction.pdf eq. (10)
void GreensTensorVacuum::integrate_k(double omega, cx_mat::fixed<3, 3> &GT,
                                     Tensor_Options fancy_complex,
                                     Weight_Options weight_function) const {
  if (fancy_complex == RE) {
    // Even though the real part of the Green's tensor is not implemented, a
    // default return value of an empty tensor was chosen, to allow for the
    // general structure of the polarizability to depend both on the real and
    // imaginary part of a given Green's tensor
    GT.zeros();
  } else if (fancy_complex == IM) {

    // Reset the tensor to store the final result
    GT.zeros();
    // Ensure that the integration limits are properly ordered
    if (omega >= 0) {

      // Numerically integrate the xx component
      auto F_xx = [=](double x) -> double {
        return this->integrand_k(x, omega, {0, 0}, fancy_complex,
                                 weight_function);
      };
      GT(0, 0) = cquad(F_xx, -omega / (1.0 + this->v), omega / (1.0 - this->v),
                       this->relerr, 0);

      // yy component
      auto F_yy = [=](double x) -> double {
        return this->integrand_k(x, omega, {1, 1}, fancy_complex,
                                 weight_function);
      };
      GT(1, 1) = cquad(F_yy, -omega / (1.0 + this->v), omega / (1.0 - this->v),
                       this->relerr, 0);

      // zz component
      GT(2, 2) = GT(1, 1);
    }
    // Switching the integration bounds for negative frequencies
    if (omega < 0) {

      // Numerically integrate the xx component
      auto F_xx = [=](double x) -> double {
        return this->integrand_k(x, omega, {0, 0}, fancy_complex,
                                 weight_function);
      };
      GT(0, 0) = -cquad(F_xx, omega / (1.0 - this->v), -omega / (1.0 + this->v),
                        this->relerr, 0);

      // yy component
      auto F_yy = [=](double x) -> double {
        return this->integrand_k(x, omega, {1, 1}, fancy_complex,
                                 weight_function);
      };
      GT(1, 1) = -cquad(F_yy, omega / (1.0 - this->v), -omega / (1.0 + this->v),
                        this->relerr, 0);

      // zz component
      GT(2, 2) = GT(1, 1);
    }
  }
}

// Implementation of the different integrands for the integration
// of the 2-d k-vector
// Ref: notes/VacuumFriction eq. (10) and (11)
double GreensTensorVacuum::integrand_k(double kv, double omega,
                                       const vec::fixed<2> &indices,
                                       Tensor_Options fancy_complex,
                                       Weight_Options weight_function) const {
  double omega_pl = (omega + kv * v);
  double omega_pl_quad = omega_pl * omega_pl;
  double xi_quad = omega_pl_quad - kv * kv;

  // Variable to store the final result
  double result = 0;

  // Only the imaginary part is implemented
  if (fancy_complex == IM) {
    // Compute the basis integrand of eq. (10)
    if (indices(0) == 0 && indices(1) == 0) {
      result = 0.5 * xi_quad;
    } else if (indices(0) == indices(1)) {
      result = 0.5 * (omega_pl_quad - xi_quad * 0.5);
    } else {
      return 0;
    }

    // Multply with the additional weight function f, the options can be found
    // in eq. (11)
    if (weight_function == KV) {
      result *= kv;
    } else if (weight_function == TEMP) {
      result /= (1.0 - exp(-beta * omega_pl));
    } else if (weight_function == KV_TEMP) {
      result *= kv / (1.0 - exp(-beta * omega_pl));
    } else if (weight_function == NON_LTE) {
      result *=
          (1. / (1. - exp(-beta * omega_pl)) - 1. / (1. - exp(-beta * omega)));
    } else if (weight_function == KV_NON_LTE) {
      result *= kv * (1. / (1. - exp(-beta * omega_pl)) -
                      1. / (1. - exp(-beta * omega)));
    }
  }

  return result;
}

double GreensTensorVacuum::omega_ch() const { return 0; }
