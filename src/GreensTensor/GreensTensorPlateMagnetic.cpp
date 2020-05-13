  //
// Created by hermasim on 20/02/2020.
//
#include "GreensTensorPlateMagnetic.h"
#include <cmath>
#include <complex>

GreensTensorPlateMagnetic::GreensTensorPlateMagnetic(
    double v, double beta, double za,
    std::shared_ptr<ReflectionCoefficients> reflection_coefficients, double delta_cut,
    vec::fixed<2> rel_err)
    : GreensTensorPlate(v, beta, za, reflection_coefficients, delta_cut,
                        rel_err) {}

GreensTensorPlateMagnetic::GreensTensorPlateMagnetic(std::string input_file)
    : GreensTensorPlate(input_file) {}

void GreensTensorPlateMagnetic::calculate_tensor(double omega, vec::fixed<2> k, cx_mat::fixed<3, 3> &GT) {
  // imaginary unit
  std::complex<double> I(0.0, 1.0);

  // load wavevectors from the struct into the corresponding variables
  double kx = k(0);
  double ky = k(1);
  double k_quad = kx * kx + ky * ky;

  // The tensor is calculated for a positive frequency omega
  // and afterwards transformed with respect to the sign of omega
  double omega_abs = std::abs(omega);
  double omega_quad = omega_abs * omega_abs;

  // kapppa is defined to have either a purely
  // positive real part or purely negatively imaginary part
  std::complex<double> kappa = sqrt(std::complex<double>(k_quad - omega_quad, 0.));
  kappa = std::complex<double>(std::abs(kappa.real()), -std::abs(kappa.imag()));

  // produce the reflection coefficients in s- and p-polarization
  std::complex<double> r_p, r_s;
  reflection_coefficients->calculate(omega_abs, kappa, r_p, r_s);

  // For an better overview and a efficient calculation, the
  // pre-factors of the p and s polarization are collected separately
  std::complex<double> pre = (2 * M_PI) * exp(-(2 * za) * kappa);
  std::complex<double> prefactor_s = pre * r_s * omega_quad / (k_quad - omega_quad);
  std::complex<double> prefactor_p = pre * r_p;

  // In the following, odd orders in ky are already omitted
  GT.zeros();
  GT(0, 0) = (prefactor_p * kx * kx + prefactor_s * ky * ky) * kappa / k_quad;
  GT(1, 1) = (prefactor_p * ky * ky + prefactor_s * kx * kx) * kappa / k_quad;
  GT(2, 2) = prefactor_p * k_quad / kappa;
  GT(2, 0) = I * prefactor_p * kx;
  GT(0, 2) = -GT(2, 0);
  GT(1, 0) =
      (-prefactor_s * ky * kx / k_quad + prefactor_p * ky * kx / k_quad) *
      kappa;
  GT(0, 1) = GT(1, 0);
  GT(2, 1) = prefactor_p * I * ky;
  GT(1, 2) = -GT(2, 1);

  // In case of negative frequencies, the tensor has to be hermitian transposed
  if (omega < 0) {
    GT = trans(GT);
  }
}

void GreensTensorPlateMagnetic::integrate_k(double omega, cx_mat::fixed<3, 3> &GT, Tensor_Options EE, Tensor_Options BE,
                                            Tensor_Options EB,
                                            Tensor_Options BB, Weight_Options weight_function) {

  // imaginary unit
  std::complex<double> I(0.0, 1.0);

  // initialize Green's tensor
  GT.zeros();

  // calculate the five non-zero elements of the Green's tensor. Here, the
  // symmetry in y direction was already applied. Thus, the integration only
  // considers twice the domain from 0 to pi.

  // the xx element
  // the xx element
  auto F_xx_R = [=](double x) -> double {
    return this->integrand_1d_k_R(x, omega, {0, 0}, EE, EB, BE, BB,
                                weight_function);
  };
  GT(0, 0) = cquad(F_xx_R, 0, 0.5 * M_PI, rel_err(1), 0) / M_PI;
  GT(0, 0) += cquad(F_xx_R, 0.5 * M_PI, M_PI, rel_err(1), 0) / M_PI;
  auto F_xx_I = [=](double x) -> double {
    return this->integrand_1d_k_I(x, omega, {0, 0}, EE, EB, BE, BB,
                                weight_function);
  };
  GT(0, 0) = I*cquad(F_xx_I, 0, 0.5 * M_PI, rel_err(1), 0) / M_PI;
  GT(0, 0) += I*cquad(F_xx_I, 0.5 * M_PI, M_PI, rel_err(1), 0) / M_PI;

  // the yy element
  auto F_yy_R = [=](double x) -> double {
    return this->integrand_1d_k_R(x, omega, {1, 1}, EE, EB, BE, BB,
                                weight_function);
  };
  GT(1, 1) = cquad(F_yy_R, 0, 0.5 * M_PI, rel_err(1), 0) / M_PI;
  GT(1, 1) += cquad(F_yy_R, 0.5 * M_PI, M_PI, rel_err(1), 0) / M_PI;
  auto F_yy_I = [=](double x) -> double {
    return this->integrand_1d_k_I(x, omega, {1, 1}, EE, EB, BE, BB,
                                weight_function);
  };
  GT(1, 1) = I*cquad(F_yy_I, 0, 0.5 * M_PI, rel_err(1), 0) / M_PI;
  GT(1, 1) += I*cquad(F_yy_I, 0.5 * M_PI, M_PI, rel_err(1), 0) / M_PI;

  // the zz element
  auto F_zz_R = [=](double x) -> double {
    return this->integrand_1d_k_R(x, omega, {2, 2}, EE, EB, BE, BB,
                                weight_function);
  };
  GT(2, 2) = cquad(F_zz_R, 0, 0.5 * M_PI, rel_err(1), 0) / M_PI;
  GT(2, 2) += cquad(F_zz_R, 0.5 * M_PI, M_PI, rel_err(1), 0) / M_PI;
  auto F_zz_I = [=](double x) -> double {
      return this->integrand_1d_k_I(x, omega, {2, 2}, EE, EB, BE, BB,
                                    weight_function);
  };
  GT(2, 2) = I*cquad(F_zz_I, 0, 0.5 * M_PI, rel_err(1), 0) / M_PI;
  GT(2, 2) += I*cquad(F_zz_I, 0.5 * M_PI, M_PI, rel_err(1), 0) / M_PI;

  // the zx element
  auto F_zx_R = [=](double x) -> double {
    return this->integrand_1d_k_R(x, omega, {2, 0}, EE, EB, BE, BB,
                              weight_function);
  };
  GT(2, 0) = cquad(F_zx_R, 0, 0.5 * M_PI, rel_err(1), 0) / M_PI;
  GT(2, 0) += cquad(F_zx_R, 0.5 * M_PI, M_PI, rel_err(1), 0) / M_PI;
  auto F_zx_I = [=](double x) -> double {
      return this->integrand_1d_k_I(x, omega, {2, 0}, EE, EB, BE, BB,
                                    weight_function);
  };
  GT(2, 0) = I*cquad(F_zx_I, 0, 0.5 * M_PI, rel_err(1), 0) / M_PI;
  GT(2, 0) += I*cquad(F_zx_I, 0.5 * M_PI, M_PI, rel_err(1), 0) / M_PI;

  // the xz
  auto F_xz_R = [=](double x) -> double {
  return this->integrand_1d_k_R(x, omega, {0, 2}, EE, EB, BE, BB,
                              weight_function);
  };
  GT(0, 2) = cquad(F_xz_R, 0, 0.5 * M_PI, rel_err(1), 0) / M_PI;
  GT(0, 2) += cquad(F_xz_R, 0.5 * M_PI, M_PI, rel_err(1), 0) / M_PI;
  auto F_xz_I = [=](double x) -> double {
      return this->integrand_1d_k_I(x, omega, {0, 2}, EE, EB, BE, BB,
                                    weight_function);
  };
  GT(0, 2) = I*cquad(F_xz_I, 0, 0.5 * M_PI, rel_err(1), 0) / M_PI;
  GT(0, 2) += I*cquad(F_xz_I, 0.5 * M_PI, M_PI, rel_err(1), 0) / M_PI;

}

double GreensTensorPlateMagnetic::integrand_1d_k_R(double phi, double omega, const vec::fixed<2> &indices,
                                                   Tensor_Options EE, Tensor_Options EB,
                                                   Tensor_Options BE, Tensor_Options BB,
                                                   Weight_Options weight_function) {
  double result = 0.;
  // import parameters
  double v = this->v;

  // The cut-off parameters acts as upper bound of the kappa integration.
  double kappa_cut = this->get_delta_cut() / (2 * this->get_za());
  // import demande relative accuracy of the integration
  vec::fixed<2> rel_err;
  rel_err(0) = this->get_rel_err_0();
  rel_err(1) = this->get_rel_err_1();
  // read integration variable phi
  double cos_phi = std::cos(phi);

  // Calculate the integrand corresponding to the given options. To resolve the
  // probably sharp edge of the Bose-Einstein distribution, the integration is
  // split at the edge, if the edged lies below the cut-off kappa_cut.
  auto F = [=](double x) -> double {
      return this->integrand_2d_k_R(x, omega, phi, indices, EE, EB, BE, BB, weight_function);
  };
  if (kappa_cut > std::abs(omega / (v * cos_phi))) {
    result = cquad(F, -std::abs(omega), 0,
                   rel_err(0), 0);
    result += cquad(F, 0,
                    std::abs(omega / (v * cos_phi)), rel_err(0), 0);
    result += cquad(F,
                    std::abs(omega / (v * cos_phi)), kappa_cut, rel_err(0), 0);
  } else {
    result = cquad(F, -std::abs(omega),
                   kappa_cut, rel_err(0), 0);
  }
  return result;
}

double GreensTensorPlateMagnetic::integrand_1d_k_I(double phi, double omega, const vec::fixed<2> &indices,
                                                   Tensor_Options EE, Tensor_Options EB,
                                                   Tensor_Options BE, Tensor_Options BB,
                                                   Weight_Options weight_function) {
  double result = 0.;
  // import parameters
  double v = this->v;

  // The cut-off parameters acts as upper bound of the kappa integration.
  double kappa_cut = this->get_delta_cut() / (2 * this->get_za());
  // import demands relative accuracy of the integration
  vec::fixed<2> rel_err;
  rel_err(0) = this->get_rel_err_0();
  rel_err(1) = this->get_rel_err_1();
  // read integration variable phi
  double cos_phi = std::cos(phi);

  // Calculate the integrand corresponding to the given options. To resolve the
  // probably sharp edge of the Bose-Einstein distribution, the integration is
  // split at the edge, if the edged lies below the cut-off kappa_cut.
  auto F = [=] (double x) -> double {
      return this->integrand_2d_k_I(x, omega, phi, indices, EE, EB, BE, BB, weight_function);
  };
  if (kappa_cut > std::abs(omega / (v * cos_phi))) {
    result = cquad(F, -std::abs(omega), 0,
                   rel_err(0), 0);
    result += cquad(F, 0,
                    std::abs(omega / (v * cos_phi)), rel_err(0), 0);
    result += cquad(F,
                    std::abs(omega / (v * cos_phi)), kappa_cut, rel_err(0), 0);
  } else {
    result = cquad(F, -std::abs(omega),
                   kappa_cut, rel_err(0), 0);
  }
  return result;
};

double GreensTensorPlateMagnetic::integrand_2d_k_R(double kappa_double, double omega, double phi,
                                                   const vec::fixed<2> &indices,
                                                   Tensor_Options EE, Tensor_Options EB, Tensor_Options BE,
                                                   Tensor_Options BB,
                                                   Weight_Options weight_function) {
  // read general input parameters
  double beta = this->get_beta();
  double v = this->get_v();
  double v_quad = v * v;
  double za = this->get_za();
  // read omega and phi from the option struct
  double omega_quad = omega * omega;
  double cos_quad = pow(cos(phi), 2);
  double sin_quad = 1.0 - cos_quad;

  // Before the real or imaginary part of the chosen matrix element can be
  // calculated, the complex result is stored in result_complex.
  std::complex<double> result_complex = (0., 0.);
  double result = 0;

  // imaginary unit
  std::complex<double> I(0.0, 1.0);

  // Transfer kappa to the correct complex value
  std::complex<double> kappa_complex = {0,0};
  double kappa_quad = 0;
  if (kappa_double < 0.0) {
    kappa_complex = std::complex<double>(0.0, kappa_double);
    kappa_quad = -kappa_double * kappa_double;
  } else {
    kappa_complex = std::complex<double>(kappa_double, 0.0);
    kappa_quad = kappa_double * kappa_double;
  }

  // Express k via frequency and kappa
  double k = (sqrt(kappa_quad * (1.0 - v_quad * cos_quad) + omega_quad) +
       v * omega * cos(phi)) /
      (1 - v_quad * cos_quad);
  double k_quad = k * k;

  // Define the Doppler-shifted frequency
  double omega_pl = (omega + k * cos(phi) * v);
  double omega_pl_quad = omega_pl * omega_pl;

  // In order to obey reality in time, a positive omega_pl is used for the
  // actual calculation. Afterwards, the corresponding symmetry operation is
  // performed if the sign of omega_pl is negative.
  double omega_pl_abs = std::abs(omega_pl);

  // producing the reflection coefficients in p- and s-polarization
  std::complex<double> r_p, r_s;
  this->reflection_coefficients->calculate(omega_pl_abs, kappa_complex, r_p, r_s);

  // General prefactor for all components
  double prefactor =
      (real(kappa_complex)-imag(kappa_complex)) / (1. - cos(phi) * v * omega_pl / k);

  // Enforce the crossing relation
  if (omega_pl < 0) {
    kappa_complex = conj(kappa_complex);
    r_s = conj(r_s);
    r_p = conj(r_p);
  }

  // prefactors for the ss-polarization, the hermitian part of the pp
  // polarization and the rotational part of the pp-polarization
  double r_s_factor;
  double r_p_herm_factor;
  std::complex<double> r_p_rot_factor;

  // check wether the electric Green's tensor should be computed
  if (EE != IGNORE) {

    if (EE == IM) {
      r_s_factor = prefactor * omega_pl_quad *
                   imag(r_s/kappa_complex * exp(-2 * za * kappa_complex));
      r_p_herm_factor =
          prefactor * imag(r_p * kappa_complex * exp(-2. * za * kappa_complex));
    }
    if (EE == RE) {
      r_s_factor = prefactor * omega_pl_quad *
                   real(r_s / kappa_complex * exp(-2 * za * kappa_complex));
      r_p_herm_factor =
          prefactor * real(r_p * kappa_complex * exp(-2. * za * kappa_complex));
    }

    // For real and complex kappa only the diagonal elements are real
      // Calculate the G_xx element
      if (indices(0) == 0 && indices(1) == 0) {
        result_complex += r_s_factor * sin_quad;
        result_complex += r_p_herm_factor * cos_quad;
      }
      // Calculate the G_yy element
      else if (indices(0) == 1 && indices(1) == 1) {
        result_complex += r_p_herm_factor * sin_quad;
        result_complex += r_s_factor * cos_quad;
      }
      // Calculate the G_zz element
      else if (indices(0) == 2 && indices(1) == 2) {
        result_complex += r_p_herm_factor * k_quad / kappa_quad;
      }
  }
  // reset the result_complex variable to store additional components of the
  // other Green's tensors check if the BE-GreensTensor should be computed
  if (BE != IGNORE) {
    if (BE == IM) {
      r_s_factor = prefactor * omega_pl_quad *
                   imag(r_s / kappa_complex * exp(-2 * za * kappa_complex));
      r_p_herm_factor =
          prefactor * imag(r_p * kappa_complex * exp(-2. * za * kappa_complex));
      r_p_rot_factor =
          prefactor * imag(r_p * exp(-2. * za * kappa_complex));
    }
    if (BE == RE) {
      r_s_factor = prefactor * omega_pl_quad *
                   real(r_s / kappa_complex * exp(-2 * za * kappa_complex));
      r_p_herm_factor =
          prefactor * real(r_p * kappa_complex * exp(-2. * za * kappa_complex));
      r_p_rot_factor =
          prefactor * real(r_p * exp(-2. * za * kappa_complex));
    }
    // Compute all the different components of the Green's tensor
    // for all kappa the diagonal elements are real
    // G_yy
    if (indices(0) == 1 && indices(1) == 1) {
      result_complex -= r_s_factor * k * v * cos(phi) / omega_pl;
    }
    // G_zz
    else if (indices(0) == 2 && indices(1) == 2) {
      result_complex -= r_p_herm_factor * pow(k, 3) * v * cos(phi) /
                        (kappa_quad * omega_pl);
    }
    //Add terms linear in kappa depending on the integration regime
    if(kappa_double >= 0) {
      if(indices(0) == 2 && indices(1) == 2) {
        result_complex += r_p_rot_factor*k*v*kappa_double*cos(phi)/omega_pl;
      }
    }
    if(kappa_double < 0)
    {
      if(indices(0) == 2 && indices(1) == 0) {
        result_complex += r_s_factor*I*v*kappa_complex*sin_quad/omega_pl;
        result_complex += r_p_herm_factor*I*v*kappa_complex*cos_quad/omega_pl;
      }
    }

  }

  // if the EB-GreensTensor should be computed
  if (EB != IGNORE) {
    // General prefactor for all components
    if (EB == IM) {
      r_s_factor = prefactor * omega_pl_quad *
                   imag(r_s / kappa_complex * exp(-2 * za * kappa_complex));
      r_p_herm_factor =
          prefactor * imag(r_p * kappa_complex * exp(-2. * za * kappa_complex));
      r_p_rot_factor =
          prefactor * imag(r_p * exp(-2. * za * kappa_complex));
    }
    if (EB == RE) {
      r_s_factor = prefactor * omega_pl_quad *
                   real(r_s / kappa_complex * exp(-2 * za * kappa_complex));
      r_p_herm_factor =
          prefactor * real(r_p * kappa_complex * exp(-2. * za * kappa_complex));
      r_p_rot_factor =
          prefactor * real(r_p * exp(-2. * za * kappa_complex));
    }
    // Compute all the different components of the Green's tensor
    // for all kappa the diagonal elements are real
    // G_yy
    if (indices(0) == 1 && indices(1) == 1) {
      result_complex -= r_s_factor * k * v * cos(phi) / omega_pl;
    }
      // G_zz
    else if (indices(0) == 2 && indices(1) == 2) {
      result_complex -= r_p_herm_factor * pow(k, 3) * v * cos(phi) /
                        (kappa_quad * omega_pl);
    }
    //Add terms linear in kappa depending on the integration regime
    if(kappa_double >= 0) {
      if(indices(0) == 2 && indices(1) == 2) {
        result_complex += r_p_rot_factor*k*v*kappa_double*cos(phi)/omega_pl;
      }
    }
    if(kappa_double < 0)
    {
      if(indices(0) == 0 && indices(1) == 2) {
        result_complex -= r_s_factor*I*v*conj(kappa_complex)*sin_quad/omega_pl;
        result_complex -= r_p_herm_factor*I*v*conj(kappa_complex)*cos_quad/omega_pl;
      }
    }

  }

if (BB != IGNORE) {
    if (BB == IM) {
      r_s_factor = prefactor * omega_pl_quad *
                   imag(r_s / kappa_complex * exp(-2 * za * kappa_complex));
      r_p_herm_factor =
          prefactor * imag(r_p * kappa_complex * exp(-2. * za * kappa_complex));
      r_p_rot_factor =
          prefactor * kappa_complex * imag(r_p * exp(-2. * za * kappa_complex));
    }
    if (BB == RE) {
      r_s_factor = prefactor * omega_pl_quad *
                   real(r_s / kappa_complex * exp(-2 * za * kappa_complex));
      r_p_herm_factor =
          prefactor * real(r_p * kappa_complex * exp(-2. * za * kappa_complex));
      r_p_rot_factor =
          prefactor * kappa_complex * real(r_p * exp(-2. * za * kappa_complex));
    }
    // Compute all the different components of the Green's tensor
    // G^BB has only real entries
    // G_yy
    if (indices(0) == 1 && indices(1) == 1) {
      result_complex += r_s_factor * k_quad * pow(v, 2) / omega_pl_quad;
    }
    // G_zz
    else if (indices(0) == 2 && indices(1) == 2) {
          result_complex +=
              r_s_factor * pow(v, 2) * kappa_quad * sin_quad / omega_pl_quad;
      if(kappa_double < 0) {
        result_complex +=
            r_p_herm_factor * pow(v, 2) *
            (pow(k, 4) + pow(kappa_complex, 3) * conj(kappa_complex)) * cos_quad /
            (kappa_quad * omega_pl_quad);
      }
      else if(kappa_double >= 0) {
        result_complex += r_p_herm_factor * pow(v,2) * omega_pl_quad * cos_quad / kappa_quad;
      }
    }
  }

  // missing I has to be added after integration procedure
  if (imag(result_complex) == 0) {
    result = real(result_complex);
  } else {
    std::cerr << "Result_complex is not purely real. There seems to be a "
                 "mistake in the implementation."
              << std::endl;
    exit(0);
  }

  // Add weighting function if demanded
  if (weight_function == KV) {
    result *= k * cos(phi);
  } else if (weight_function == TEMP) {
    result /= (1.0 - exp(-beta * omega_pl));
  } else if (weight_function == NON_LTE) {
    result *=
        1. / (1.0 - exp(-beta * omega_pl)) - 1. / (1.0 - exp(-beta * omega));
  } else if (weight_function == KV_TEMP) {
    result *= k * cos(phi) / (1.0 - exp(-beta * omega_pl));
  } else if (weight_function == KV_NON_LTE) {
    result *=
        k * cos(phi) *
        (1. / (1.0 - exp(-beta * omega_pl)) - 1. / (1.0 - exp(-beta * omega)));
  }

  return result;
}

double GreensTensorPlateMagnetic::integrand_2d_k_I(double kappa_double, double omega, double phi,
                                                   const vec::fixed<2> &indices,
                                                   Tensor_Options EE, Tensor_Options EB, Tensor_Options BE,
                                                   Tensor_Options BB,
                                                   Weight_Options weight_function) {
  // read general input parameters
  double beta = this->beta;
  double v = this->v;
  double v_quad = v * v;
  double za = this->za;
  // read omega and phi from the option struct
  double omega_quad = omega * omega;
  double cos_quad = pow(cos(phi), 2);
  double sin_quad = 1.0 - cos_quad;

  // Before the real or imaginary part of the chosen matrix element can be
  // calculated, the complex result is stored in result_complex.
  std::complex<double> result_complex = {0., 0.};
  double result = 0;

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

  // Express k via frequency and kappa
  double k = (sqrt(kappa_quad * (1.0 - v_quad * cos_quad) + omega_quad) +
       v * omega * cos(phi)) /
      (1 - v_quad * cos_quad);
  double k_quad = k * k;

  // Define the Doppler-shifted frequency
  double omega_pl = (omega + k * cos(phi) * v);
  double omega_pl_quad = omega_pl * omega_pl;

  // In order to obey reality in time, a positive omega_pl is used for the
  // actual calculation. Afterwards, the corresponding symmetry operation is
  // performed if the sign of omega_pl is negative.
  double omega_pl_abs = std::abs(omega_pl);

  // reflection coefficients and pre-factors of the corresponding polarization
  std::complex<double> r_p, r_s;
  // producing the reflection coefficients in p- and s-polarization
  reflection_coefficients->ref(r_p, r_s, omega_pl_abs, kappa_complex);

  //General prefactor for all components
  double prefactor =
      (real(kappa_complex)-imag(kappa_complex)) / (1. - cos(phi) * v * omega_pl / k);

  // Enforce the crossing relation
  if (omega_pl < 0) {
    kappa_complex = conj(kappa_complex);
    r_s = conj(r_s);
    r_p = conj(r_p);
  }


  // prefactors for the ss-polarization, the hermitian part of the pp
  // polarization and the rotational part of the pp-polarization
  double r_s_factor;
  double r_p_herm_factor;
  std::complex<double> r_p_rot_factor;

  // check wether the electric Green's tensor should be computed
  if (EE != IGNORE) {
    if (EE == IM) {
      r_p_rot_factor =
          prefactor * imag(r_p * exp(-2. * za * kappa_complex));
    }
    if (EE == RE) {
      r_p_rot_factor =
          prefactor * real(r_p * exp(-2. * za * kappa_complex));
    }

    // For real and complex kappa only the off-diagonal elements are complex
      // Calculate the G_zx element
      if (indices(0) == 2 && indices(1) == 0) {
        result_complex += r_p_rot_factor * I * k * cos(phi) ;
      }
      // Calculate the G_xz element
      else if (indices(0) == 0 && indices(1) == 2) {
        result_complex += -r_p_rot_factor * I * k * cos(phi);
      }
  }
  // check if the BE-GreensTensor should be computed
  if (BE != IGNORE) {
    if (BE == IM) {
      r_s_factor = prefactor * omega_pl_quad *
                   imag(r_s / kappa_complex * exp(-2 * za * kappa_complex));
      r_p_herm_factor =
          prefactor * imag(r_p * kappa_complex * exp(-2. * za * kappa_complex));
      r_p_rot_factor =
          prefactor * imag(r_p * exp(-2. * za * kappa_complex));
    }
    if (BE == RE) {
      r_s_factor = prefactor * omega_pl_quad *
                   real(r_s / kappa_complex * exp(-2 * za * kappa_complex));
      r_p_herm_factor =
          prefactor * real(r_p * kappa_complex * exp(-2. * za * kappa_complex));
      r_p_rot_factor =
          prefactor * real(r_p * exp(-2. * za * kappa_complex));
    }
    // Compute all the different components of the Green's tensor
    // for real kappa the diagonal elements are complex
      // G_zx
    if (indices(0) == 2 && indices(1) == 0) {
      result_complex -= r_p_rot_factor * I * k_quad * v * cos_quad /
                        omega_pl;
    }
    //Add terms depending on the integration regime
    if(kappa_double >= 0) {
      if (indices(0) == 2 && indices(1) == 0) {
        result_complex +=
            r_s_factor * I * v * kappa_complex * sin_quad / omega_pl;
        result_complex +=
            r_p_herm_factor * I * v * kappa_complex * cos_quad / omega_pl;
      }
    }
    if(kappa_double < 0) {
      if(indices(0) == 2 && indices(1) == 2) {
        result_complex += r_p_rot_factor*k*v*kappa_complex*cos(phi)/omega_pl;
      }
    }
  }
  //Check if G^BE should be computed
  if (EB != IGNORE) {
    if (EB == IM) {
      r_s_factor = prefactor * omega_pl_quad *
                   imag(r_s / kappa_complex * exp(-2 * za * kappa_complex));
      r_p_herm_factor =
          prefactor * imag(r_p * kappa_complex * exp(-2. * za * kappa_complex));
      r_p_rot_factor =
          prefactor * imag(r_p * exp(-2. * za * kappa_complex));
    }
    if (EB == RE) {
      r_s_factor = prefactor * omega_pl_quad *
                   real(r_s / kappa_complex * exp(-2 * za * kappa_complex));
      r_p_herm_factor =
          prefactor * real(r_p * kappa_complex * exp(-2. * za * kappa_complex));
      r_p_rot_factor =
          prefactor * real(r_p * exp(-2. * za * kappa_complex));
    }
    // Compute all the different components of the Green's tensor
    // for real kappa the diagonal elements are complex
    // G_zx
    if (indices(0) == 0 && indices(1) == 2) {
      result_complex += r_p_rot_factor * I * k_quad * v * cos_quad /
                        omega_pl;
    }

    //Add terms depending on the integration regime
    if(kappa_double >= 0) {
      if (indices(0) == 0 && indices(1) == 2) {
        result_complex -=
            r_s_factor * I * v * kappa_double * sin_quad / omega_pl;
        result_complex -=
            r_p_herm_factor * I * v * kappa_double * cos_quad / omega_pl;
      }
    }
    if(kappa_double < 0) {
      if(indices(0) == 2 && indices(1) == 2) {
        result_complex -= r_p_rot_factor*k*v*kappa_complex*cos(phi)/omega_pl;
      }
    }
  }

  // There are no complex contribution of the G^BB Green"s tensor.

  // Ensure that the integrand was implemented correctly
  if (real(result_complex) == 0) {
    result = imag(result_complex);
  } else {
    std::cerr << "Result_complex is not purely imaginary. There seems to be a "
                 "mistake in the implementation."
              << std::endl;
    exit(0);
  }

  // Add weighting function if demanded
  if (weight_function == KV) {
    result *= k * cos(phi);
  } else if (weight_function == TEMP) {
    result /= (1.0 - exp(-beta * omega_pl));
  } else if (weight_function == NON_LTE) {
    result *=
        1. / (1.0 - exp(-beta * omega_pl)) - 1. / (1.0 - exp(-beta * omega));
  } else if (weight_function == KV_TEMP) {
    result *= k * cos(phi) / (1.0 - exp(-beta * omega_pl));
  } else if (weight_function == KV_NON_LTE) {
    result *=
        k * cos(phi) *
        (1. / (1.0 - exp(-beta * omega_pl)) - 1. / (1.0 - exp(-beta * omega)));
  }
  return result;
}
