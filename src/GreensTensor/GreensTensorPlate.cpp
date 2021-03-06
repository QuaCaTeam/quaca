// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <utility>
namespace pt = boost::property_tree;

// integration routine
#include "../Calculations/Integrations.h"

#include "../Permittivity/PermittivityFactory.h"
#include "../ReflectionCoefficients/ReflectionCoefficientsFactory.h"
#include "GreensTensorPlate.h"

GreensTensorPlate::GreensTensorPlate(
    double v, double beta, double za,
    std::shared_ptr<ReflectionCoefficients> reflection_coefficients,
    double delta_cut, const vec::fixed<2> &rel_err)
    : GreensTensor(v, beta), za(za), delta_cut(delta_cut), rel_err(rel_err),
      reflection_coefficients(std::move(reflection_coefficients)) {

  // assertions
  assert(this->za >= 0);
  assert(this->delta_cut >= 0);
  assert(this->rel_err(0) >= 0 && this->rel_err(1) >= 0);
}

GreensTensorPlate::GreensTensorPlate(const std::string &input_file)
    : GreensTensor(input_file) {
  this->reflection_coefficients =
      ReflectionCoefficientsFactory::create(input_file);

  // Create a root
  pt::ptree root;

  // Load the json file in this ptree
  pt::read_json(input_file, root);

  // check if type is right
  std::string type = root.get<std::string>("GreensTensor.type");
  assert(type == "plate");

  // read parameters
  this->za = root.get<double>("GreensTensor.za");
  this->delta_cut = root.get<double>("GreensTensor.delta_cut");
  this->rel_err(0) = root.get<double>("GreensTensor.rel_err_0");
  this->rel_err(1) = root.get<double>("GreensTensor.rel_err_1");

  // assertions
  assert(this->za >= 0);
  assert(this->delta_cut >= 0);
  assert(this->rel_err(0) >= 0 && this->rel_err(1) >= 0);
}

void GreensTensorPlate::calculate_tensor(double omega, vec::fixed<2> k,
                                         cx_mat::fixed<3, 3> &GT) const {
  // imaginary unit
  std::complex<double> I(0.0, 1.0);

  // load wavevectors from the struct into the corresponding variables
  double kx = k(0);
  double ky = k(1);
  double k_quad = kx * kx + ky * ky;

  // squared and absolute value of the frequency
  // The tensor is calculated for a positive frequency omega
  // and afterwards transformed with respect to the sign of omega
  double omega_abs = std::abs(omega);
  double omega_quad = omega_abs * omega_abs;

  // propagation in free space and within the surface material
  // kappa is defined to have either a purely
  // positive real part or purely negatively imaginary part
  std::complex<double> kappa =
      sqrt(std::complex<double>(k_quad - omega_quad, 0.));
  kappa = std::complex<double>(std::abs(kappa.real()), -std::abs(kappa.imag()));

  // produce the reflection coefficients in s- and p-polarization
  std::complex<double> r_p, r_s;
  reflection_coefficients->calculate(omega_abs, kappa, r_p, r_s);

  // For an better overview and a efficient calculation, the
  // pre-factors of the p and s polarization are collected separately
  std::complex<double> pre = (2 * M_PI) * exp(-(2 * za) * kappa);
  std::complex<double> prefactor_s =
      pre * r_s * omega_quad / (k_quad - omega_quad);
  std::complex<double> prefactor_p = pre * r_p;

  // In the following, odd orders in ky are already omitted
  GT.zeros();
  GT(0, 0) = (prefactor_p * kx * kx + prefactor_s * ky * ky) * kappa / k_quad;
  GT(1, 1) = (prefactor_p * ky * ky + prefactor_s * kx * kx) * kappa / k_quad;
  GT(2, 2) = prefactor_p * k_quad / kappa;
  GT(2, 0) = I * prefactor_p * kx;
  GT(0, 2) = -GT(2, 0);

  // In case of negative frequencies, the tensor has to be hermitian transposed
  if (omega < 0) {
    GT = trans(GT);
  }
}

void GreensTensorPlate::integrate_k(double omega, cx_mat::fixed<3, 3> &GT,
                                    Tensor_Options fancy_complex,
                                    Weight_Options weight_function) const {

  // imaginary unit
  std::complex<double> I(0.0, 1.0);

  // intialize Green's tensor
  GT.zeros();

  // calculate the five non-zero elements of the Green's tensor. Here, the
  // symmetry in y direction was already applied. Thus, the integration only
  // consideres twice the domain from 0 to pi.

  // the xx element
  auto F_xx = [=](double x) -> double {
    return this->integrand_1d_k(x, omega, {0, 0}, fancy_complex,
                                weight_function);
  };
  GT(0, 0) = cquad(F_xx, 0, M_PI, rel_err(1), 0) / M_PI;

  // the yy element
  auto F_yy = [=](double x) -> double {
    return this->integrand_1d_k(x, omega, {1, 1}, fancy_complex,
                                weight_function);
  };
  GT(1, 1) = cquad(F_yy, 0, M_PI, rel_err(1), 0) / M_PI;

  // the zz element
  auto F_zz = [=](double x) -> double {
    return this->integrand_1d_k(x, omega, {2, 2}, fancy_complex,
                                weight_function);
  };
  GT(2, 2) = cquad(F_zz, 0, M_PI, rel_err(1), 0) / M_PI;

  // the zx element
  auto F_zx = [=](double x) -> double {
    return this->integrand_1d_k(x, omega, {2, 0}, fancy_complex,
                                weight_function);
  };
  GT(2, 0) = I * cquad(F_zx, 0, M_PI, rel_err(1), 0) / M_PI;

  // the xz element
  GT(0, 2) = -GT(2, 0);
}

double GreensTensorPlate::integrand_1d_k(double phi, double omega,
                                         const uvec::fixed<2> &indices,
                                         Tensor_Options fancy_complex,
                                         Weight_Options weight_function) const {

  double result;

  // The cut-off parameters acts as upper bound of the kappa integration.
  double kappa_cut = delta_cut / (2 * za);

  // read integration variable phi
  double cos_phi = std::cos(phi);

  // define integrand
  auto F = [=](double x) -> double {
    return this->integrand_2d_k(x, omega, phi, indices, fancy_complex,
                                weight_function);
  };

  // Calculate low-temperature edge
  double edge = std::abs(omega / (v * cos_phi));

  // Calculate the integrand corresponding to the given options. To resolve the
  // probably sharp edge of the Bose-Einstein distribution, the integration is
  // split at the edge, if the edged lies below the cut-off kappa_cut.
  if ( (kappa_cut > edge) && (2*za/v < beta) ) {
    result = cquad(F, 0, std::abs(omega / (v * cos_phi)), rel_err(0), 0);
    result += cquad(F, edge, kappa_cut, rel_err(0), std::abs(result)*rel_err(0));
  } else {
    result = cquad(F, 0, kappa_cut, rel_err(0), 0);
  }
    result += cquad(F, -std::abs(omega), 0, rel_err(0), std::abs(result)*rel_err(0));

  return result;
}

double GreensTensorPlate::integrand_2d_k(double kappa_double, double omega,
                                         double phi,
                                         const uvec::fixed<2> &indices,
                                         Tensor_Options fancy_complex,
                                         Weight_Options weight_function) const {

  double v_quad = v * v;
  double omega_quad = omega * omega;
  double cos_phi = cos(phi);
  double cos_phi_quad = cos_phi * cos_phi;
  double sin_phi_quad = 1.0 - cos_phi_quad;

  // Before the real or imaginary part of the chosen matrix element can be
  // calculated, the complex result is stored in result_complex.
  std::complex<double> result_complex;
  double result = 0.;

  // imaginary unit
  std::complex<double> I(0.0, 1.0);

  // permittivity and propagation through vacuum (kappa) and surface material
  std::complex<double> kappa_complex;
  double kappa_quad;
  // Transfer kappa to the correct complex value
  if (kappa_double < 0.0) {
    kappa_complex = std::complex<double>(0.0, kappa_double);
    kappa_quad = -kappa_double * kappa_double;
  } else {
    kappa_complex = std::complex<double>(kappa_double, 0.0);
    kappa_quad = kappa_double * kappa_double;
  }

  // Express kappa via frequency and kappa.
  // In order to achieve the desired accuracy, we subtract first 
  // (kappa^2 + omega^2), since this might be equal zero.
  double k = (sqrt((kappa_quad + omega_quad)- kappa_quad * v_quad * cos_phi_quad) +
              v * omega * cos_phi) / (1.E0 - v_quad * cos_phi_quad);
  double k_quad = k * k;

  // Define the Doppler-shifted frequency
  double omega_pl = (omega + k * cos_phi * v);
  double omega_pl_quad = omega_pl * omega_pl;

  // In order to obey reality in time, a positive omega_pl is used for the
  // actual calculation. Afterwards, the corresponding symmetry operation is
  // performed if the sign of omega_pl is negative.
  double omega_pl_abs = std::abs(omega_pl);

  // producing the reflection coefficients in p- and s-polarization
  // reflection coefficients and pre-factors of the corresponding polarization
  std::complex<double> r_p, r_s;
  reflection_coefficients->calculate(omega_pl_abs, kappa_complex, r_p, r_s);

  // Impose reality in time
  if (omega_pl < 0) {
    r_s = conj(r_s);
    r_p = conj(r_p);
    kappa_complex = conj(kappa_complex);
  }

  // helpful prefactors
  // general prefactor with volume element and exponential
  std::complex<double> prefactor = std::abs(kappa_complex) *
                                   exp(-2 * za * kappa_complex) /
                                   (1. - cos_phi * v * omega_pl / k);
  // For an better overview and a efficient calculation, we collect the
  // pre-factors of the p and s polarization separately
  std::complex<double> prefactor_s =
      prefactor * r_s * omega_pl_quad / kappa_complex;
  std::complex<double> prefactor_p = prefactor * r_p * kappa_complex;

  // Calculate the G_xx element
  if (indices(0) == 0 && indices(1) == 0) {
    result_complex = prefactor_p * cos_phi_quad + prefactor_s * sin_phi_quad;
  }
  // Calculate the G_yy element
  else if (indices(0) == 1 && indices(1) == 1) {
    result_complex = prefactor_p * sin_phi_quad + prefactor_s * cos_phi_quad;
  }
  // Calculate the G_zz element
  else if (indices(0) == 2 && indices(1) == 2) {
    result_complex = prefactor_p * k_quad / kappa_quad;
  }
  // Calculate the G_zx element
  else if (indices(0) == 2 && indices(1) == 0) {
    result_complex = prefactor_p * I * cos_phi * k / kappa_complex;
  }
  // Calculate the G_xz element
  else if (indices(0) == 0 && indices(1) == 2) {
    result_complex = -prefactor_p * I * cos_phi * k / kappa_complex;
  } else {
    result_complex = 0.;
  }

  // Add weighting function if demanded
  if (weight_function == KV) {
    result_complex *= k * cos_phi;
  } else if (weight_function == TEMP) {
    result_complex /= (1.0 - exp(-beta * omega_pl));
  } else if (weight_function == NON_LTE) {
    result_complex *=
        1. / (1.0 - exp(-beta * omega_pl)) - 1. / (1.0 - exp(-beta * omega));
  } else if (weight_function == KV_TEMP) {
    result_complex *= k * cos_phi / (1.0 - exp(-beta * omega_pl));
  } else if (weight_function == KV_NON_LTE) {
    result_complex *=
        k * cos_phi *
        (1. / (1.0 - exp(-beta * omega_pl)) - 1. / (1.0 - exp(-beta * omega)));
  }

  // Calculate fancy real part of the given matrix element
  if (fancy_complex == RE) {
    if ((indices(0) == 2 && indices(1) == 0) ||
        (indices(0) == 0 && indices(1) == 2)) {
      // Mind the missing leading I! This must be added after the double
      // integration!
      result = result_complex.imag();
    } else {
      result = result_complex.real();
    }
  }
  // Calculate fancy imaginary part of the given matrix element
  else if (fancy_complex == IM) {
    if ((indices(0) == 2 && indices(1) == 0) ||
        (indices(0) == 0 && indices(1) == 2)) {
      // Mind the missing leading I! This must be added after the double
      // integration!
      result = -result_complex.real();
    } else {
      result = result_complex.imag();
    }
  }

  return result;
}

std::complex<double> GreensTensorPlate::get_r_p(double omega, double k) const {
  std::complex<double> r_p, r_s;
  std::complex<double> kappa;
  if (k < omega) {
    kappa = std::complex<double>(0., -sqrt(omega * omega - k * k));
  } else {
    kappa = std::complex<double>(sqrt(k * k - omega * omega), 0.);
  }
  reflection_coefficients->calculate(omega, kappa, r_p, r_s);
  return r_p;
}

std::complex<double> GreensTensorPlate::get_r_s(double omega, double k) const {
  std::complex<double> r_s, r_p;
  std::complex<double> kappa;
  if (k < omega) {
    kappa = std::complex<double>(0., -sqrt(omega * omega - k * k));
  } else {
    kappa = std::complex<double>(sqrt(k * k - omega * omega), 0.);
  }
  reflection_coefficients->calculate(omega, kappa, r_p, r_s);
  return r_s;
}

double GreensTensorPlate::omega_ch() const {
  // Calculate omega_cut (reasonable for every plate setup)
  return this->delta_cut * this->v / this->za;
}

void GreensTensorPlate::print_info(std::ostream &stream) const {
  stream << "# GreensTensorPlate\n#\n"
         << "# v = " << v << "\n"
         << "# beta = " << beta << "\n"
         << "# za = " << za << "\n"
         << "# delta_cut = " << delta_cut << "\n"
         << "# rel_err = " << rel_err(0) << "," << rel_err(1) << "\n";
}
