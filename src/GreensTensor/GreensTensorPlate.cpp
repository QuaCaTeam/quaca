// json parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

// integration routine
#include "../Calculations/Integrations.h"
#include "../Permittivity/PermittivityFactory.h"
#include "../ReflectionCoefficients/ReflectionCoefficientsFactory.h"
#include "GreensTensorPlate.h"

GreensTensorPlate::GreensTensorPlate(
    double v, double za, double beta,
    ReflectionCoefficients *reflection_coefficients, double delta_cut,
    vec::fixed<2> rel_err)
    : GreensTensor(v, beta) {
  this->za = za;
  this->delta_cut = delta_cut;
  this->rel_err = rel_err;
  this->reflection_coefficients = reflection_coefficients;
}

GreensTensorPlate::GreensTensorPlate(std::string input_file)
    : GreensTensor(input_file) {
  this->reflection_coefficients =
      ReflectionCoefficientsFactory::create(input_file);

  // Create a root
  pt::ptree root;

  // Load the json file in this ptree
  pt::read_json(input_file, root);

  // check if type is right
  auto type = root.get<std::string>("GreensTensor.type");
  assert(type == "plate");

  // read parameters
  this->za = root.get<double>("GreensTensor.za");
  this->delta_cut = root.get<double>("GreensTensor.delta_cut");
  this->rel_err(0) = root.get<double>("GreensTensor.rel_err_0");
  this->rel_err(1) = root.get<double>("GreensTensor.rel_err_1");
}

void GreensTensorPlate::calculate_tensor(cx_mat::fixed<3, 3> &GT,
                                         Options_GreensTensor opts) {
  // wavevectors
  double kx, ky, k_quad;
  // squared and absolute value of the frequency
  double omega_quad, omega_abs;
  // imaginary unit
  std::complex<double> I(0.0, 1.0);
  // propagation in free space and within the surface material
  std::complex<double> kappa;
  // polarization dependent reflection coefficients
  std::complex<double> r_p, r_s;
  // some prefactors for a better overview
  std::complex<double> pre, prefactor_p, prefactor_s;

  // load wavevectors from the struct into the corresponding variables
  kx = opts.kvec(0);
  ky = opts.kvec(1);
  k_quad = kx * kx + ky * ky;

  // The tensor is calculated for a positive frequency omega
  // and afterwards transformed with respect to the sign of omega
  omega_abs = std::abs(opts.omega);
  omega_quad = omega_abs * omega_abs;

  // kapppa is defined to have either a purely
  // positive real part or purely negatively imaginary part
  kappa = sqrt(std::complex<double>(k_quad - omega_quad, 0.));
  kappa = std::complex<double>(std::abs(kappa.real()), -std::abs(kappa.imag()));

  // produce the reflection coefficients in s- and p-polarization
  reflection_coefficients->ref(r_p, r_s, omega_abs, kappa);

  // For an better overview and a efficient calculation, the
  // pre-factors of the p and s polarization are collected separately
  pre = (2 * M_PI) * exp(-(2 * za) * kappa);
  prefactor_s = pre * r_s * omega_quad / (k_quad - omega_quad);
  prefactor_p = pre * r_p;

  // In the following, odd orders in ky are already omitted
  GT.zeros();
  GT(0, 0) = (prefactor_p * kx * kx + prefactor_s * ky * ky) * kappa / k_quad;
  GT(1, 1) = (prefactor_p * ky * ky + prefactor_s * kx * kx) * kappa / k_quad;
  GT(2, 2) = prefactor_p * k_quad / kappa;
  GT(2, 0) = I * prefactor_p * kx;
  GT(0, 2) = -GT(2, 0);

  // In case of negative frequencies, the tensor has to be hermitian transposed
  if (opts.omega < 0) {
    GT = trans(GT);
  }
}

void GreensTensorPlate::integrate_k(cx_mat::fixed<3, 3> &GT,
                                    Options_GreensTensor opts) {

  // imaginary unit
  std::complex<double> I(0.0, 1.0);
  // intialize Green's tensor
  GT.zeros();

  // calculate the five non-zero elements of the Green's tensor. Here, the
  // symmetry in y direction was already applied. Thus, the integration only
  // consideres twice the domain from 0 to pi.

  // the xx element
  opts.indices = {0, 0};
  GT(0, 0) =
      cquad(&integrand_1d_k, &opts, 0, 0.5 * M_PI, this->rel_err(1), 0) / M_PI;
  GT(0, 0) +=
      cquad(&integrand_1d_k, &opts, 0.5 * M_PI, M_PI, this->rel_err(1), 0) /
      M_PI;
  // the yy element
  opts.indices = {1, 1};
  GT(1, 1) =
      cquad(&integrand_1d_k, &opts, 0, 0.5 * M_PI, this->rel_err(1), 0) / M_PI;
  GT(1, 1) +=
      cquad(&integrand_1d_k, &opts, 0.5 * M_PI, M_PI, this->rel_err(1), 0) /
      M_PI;
  // the zz element
  opts.indices = {2, 2};
  GT(2, 2) =
      cquad(&integrand_1d_k, &opts, 0, 0.5 * M_PI, this->rel_err(1), 0) / M_PI;
  GT(2, 2) +=
      cquad(&integrand_1d_k, &opts, 0.5 * M_PI, M_PI, this->rel_err(1), 0) /
      M_PI;
  // the zx element
  opts.indices = {2, 0};
  GT(2, 0) = I *
             cquad(&integrand_1d_k, &opts, 0, 0.5 * M_PI, this->rel_err(1), 0) /
             M_PI;
  GT(2, 0) +=
      I * cquad(&integrand_1d_k, &opts, 0.5 * M_PI, M_PI, this->rel_err(1), 0) /
      M_PI;
  // the xz element
  GT(0, 2) = -GT(2, 0);
}

double GreensTensorPlate::integrand_1d_k(double phi, void *opts) {
  // The needed parameters for the integration are encoded in the void pointer.
  // This void pointer is casted the options struct given in GreensTensor.h.
  auto *opts_pt = static_cast<Options_GreensTensor *>(opts);
  // To access attributes of the GreensTensorPlate class, the class pointer
  // within the struct is casted.
  auto *pt = dynamic_cast<GreensTensorPlate *>(opts_pt->class_pt);

  double result = 0;
  // import parameters
  double omega = opts_pt->omega;
  double v = pt->v;
  // The cut-off parameters acts as upper bound of the kappa integration.
  double kappa_cut = pt->delta_cut / (2 * pt->za);
  // read integration variable phi
  double cos_phi = std::cos(phi);

  // Write the integration variable into the options struct
  opts_pt->kvec(0) = phi;

  // Calculate the integrand corresponding to the given options. To resolve the
  // probably sharp edge of the Bose-Einstein distribution, the integration is
  // split at the edge, if the edged lies below the cut-off kappa_cut.
  if (kappa_cut > std::abs(omega / (v * cos_phi))) {
    //result =
     //   cquad(&integrand_2d_k, opts, -std::abs(omega), 0, pt->rel_err(0), 0);
    result += cquad(&integrand_2d_k, opts, 0, std::abs(omega / (v * cos_phi)),
                    pt->rel_err(0), 0);
    result += cquad(&integrand_2d_k, opts, std::abs(omega / (v * cos_phi)),
                    kappa_cut, pt->rel_err(0), 0);
  } else {
    result = cquad(&integrand_2d_k, opts, -std::abs(omega), kappa_cut,
                   pt->rel_err(0), 0);
  }
  return result;
}

double GreensTensorPlate::integrand_2d_k(double kappa_double, void *opts) {
  // The needed parameters for the integration are encoded in the void pointer.
  // This void pointer is casted the options struct given in GreensTensor.h.
  auto *opts_pt = static_cast<Options_GreensTensor *>(opts);
  // To access attributes of the GreensTensorPlate class, the class pointer
  // within the struct is casted.
  auto *pt = dynamic_cast<GreensTensorPlate *>(opts_pt->class_pt);

  // read general input parameters
  double beta = pt->beta;
  double v = pt->v;
  double v_quad = v * v;
  double za = pt->za;
  // read omega and phi from the option struct
  double omega = opts_pt->omega;
  double omega_quad = omega * omega;
  double cos_phi = cos(opts_pt->kvec(0));
  double cos_phi_quad = cos_phi * cos_phi;
  double sin_phi_quad = 1.0 - cos_phi_quad;

  // Before the real or imaginary part of the chosen matrix element can be
  // calculated, the complex result is stored in result_complex.
  std::complex<double> result_complex;
  double result = 0;
  // Doppler-shifted frequency
  double omega_pl, omega_pl_quad;
  // wavevector
  double k, k_quad;

  // imaginary unit
  std::complex<double> I(0.0, 1.0);

  // permittivity and propagation through vacuum (kappa) and surface material
  std::complex<double> kappa_complex;
  double kappa_quad;

  // reflection coefficients and pre-factors of the corresponding polarization
  std::complex<double> r_p, r_s;

  // helpful prefactors
  std::complex<double> prefactor_p, prefactor_s, prefactor_off;

  // Transfer kappa to the correct complex value
  if (kappa_double < 0.0) {
    kappa_complex = std::complex<double>(0.0, kappa_double);
    kappa_quad = -kappa_double * kappa_double;
  } else {
    kappa_complex = std::complex<double>(kappa_double, 0.0);
    kappa_quad = kappa_double * kappa_double;
  }

  // Express k via frequency and kappa
  k = (sqrt(kappa_quad * (1.0 - v_quad * cos_phi_quad) + omega_quad) +
       v * omega * cos_phi) /
      (1 - v_quad * cos_phi_quad);
  k_quad = k * k;

  // Define the Doppler-shifted frequency
  omega_pl = (omega + k * cos_phi * v);
  omega_pl_quad = omega_pl * omega_pl;

  // In order to obey reality in time, a positive omega_pl is used for the
  // actual calculation. Afterwards, the corresponding symmetry operation is
  // performed if the sign of omega_pl is negative.
  double omega_pl_abs = std::abs(omega_pl);

  // producing the reflection coefficients in p- and s-polarization
  pt->reflection_coefficients->ref(r_p, r_s, omega_pl_abs, kappa_complex);
  // For an better overview and a efficient calculation, we collect the
  // pre-factors of the p and s polarization separately
  prefactor_p =
      exp(-2 * za * kappa_complex) / (1. - cos_phi * v * omega_pl / k);
  prefactor_s = prefactor_p * r_s * omega_pl_quad;
  prefactor_p = prefactor_p * r_p * kappa_quad;
  prefactor_off = prefactor_p * k * cos_phi / kappa_complex;

  // Impose reality in time
  if (omega_pl < 0) {
    prefactor_s = conj(prefactor_s);
    prefactor_p = conj(prefactor_p);
    prefactor_off = conj(prefactor_off);
  }

  // Calculate the G_xx element
  if (opts_pt->indices(0) == 0 && opts_pt->indices(1) == 0) {
    result_complex = prefactor_p * cos_phi_quad + prefactor_s * sin_phi_quad;
  }
  // Calculate the G_yy element
  else if (opts_pt->indices(0) == 1 && opts_pt->indices(1) == 1) {
    result_complex = prefactor_p * sin_phi_quad + prefactor_s * cos_phi_quad;
  }
  // Calculate the G_zz element
  else if (opts_pt->indices(0) == 2 && opts_pt->indices(1) == 2) {
    result_complex = prefactor_p * k_quad / kappa_quad;
  }
  // Calculate the G_zx element
  else if (opts_pt->indices(0) == 2 && opts_pt->indices(1) == 0) {
    result_complex = I * prefactor_off;
  }
  // Calculate the G_xz element
  else if (opts_pt->indices(0) == 0 && opts_pt->indices(1) == 2) {
    result_complex = -I * prefactor_off;
  } else {
    result_complex = (0, 0);
  }

  // Add weighting function if demanded
  if (opts_pt->weight_function == KV) {
    result_complex *= k * cos_phi;
  } else if (opts_pt->weight_function == TEMP) {
    result_complex /= (1.0 - exp(-beta * omega_pl));
  } else if (opts_pt->weight_function == NON_LTE) {
    result_complex *=
        1. / (1.0 - exp(-beta * omega_pl)) - 1. / (1.0 - exp(-beta * omega));
  } else if (opts_pt->weight_function == KV_TEMP) {
    result_complex *= k * cos_phi / (1.0 - exp(-beta * omega_pl));
  } else if (opts_pt->weight_function == KV_NON_LTE) {
    result_complex *=
        k * cos_phi *
        (1. / (1.0 - exp(-beta * omega_pl)) - 1. / (1.0 - exp(-beta * omega)));
  }

  // Calculate fancy real part of the given matrix element
  if (opts_pt->fancy_complex == RE) {
    if (opts_pt->indices(0) == 2 && opts_pt->indices(1) == 0 ||
        opts_pt->indices(0) == 0 && opts_pt->indices(1) == 2) {
      // Mind the missing leading I! This must be added after the double
      // integration!
      result = result_complex.imag();
    } else {
      result = result_complex.real();
    }
  }
  // Calculate fancy imaginary part of the given matrix element
  else if (opts_pt->fancy_complex == IM) {
    if (opts_pt->indices(0) == 2 && opts_pt->indices(1) == 0 ||
        opts_pt->indices(0) == 0 && opts_pt->indices(1) == 2) {
      // Mind the missing leading I! This must be added after the double
      // integration!
      result = -result_complex.real();
    } else {
      result = result_complex.imag();
    }
  }
  return result;
}

std::complex<double> GreensTensorPlate::get_r_p(double omega, double k) {
  std::complex<double> r_p, r_s;
  std::complex<double> kappa;
  if (omega > 0) {
    kappa = std::complex<double>(0., -sqrt(omega * omega - k * k));
  } else {
    kappa = std::complex<double>(sqrt(k * k - omega * omega), 0.);
  }
  reflection_coefficients->ref(r_p, r_s, omega, kappa);
  return r_p;
}

std::complex<double> GreensTensorPlate::get_r_s(double omega, double k) {
  std::complex<double> r_s, r_p;
  std::complex<double> kappa;
  if (omega > 0) {
    kappa = std::complex<double>(0., -sqrt(omega * omega - k * k));
  } else {
    kappa = std::complex<double>(sqrt(k * k - omega * omega), 0.);
  }
  reflection_coefficients->ref(r_p, r_s, omega, kappa);
  return r_s;
}

double GreensTensorPlate::omega_ch() {
  // Calculate omega_cut (reasonable for every plate setup)
  return this->delta_cut * this->v / this->za;
}
