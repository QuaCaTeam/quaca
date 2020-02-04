#include "Quaca.h"
#include "catch.hpp"
#include <armadillo>
#include <complex>

TEST_CASE("Constructors of the reflection coefficients work properly", "[ReflectionCoefficients]"){

  SECTION("Constructor with argument list works"){
    double omega_p = 9;
    double omega = 1;
    double gamma = 3.5E-2;
      Permittivity *perm = new PermittivityDrude(gamma, omega_p);
    ReflectionCoefficientsLocBulk RefC(perm);

    REQUIRE(Approx(RefC.get_epsilon(omega).real()).epsilon(1E-6) == perm->epsilon(omega).real());
  }

  SECTION("Constructor with ini file works") {
    double omega = 1;
    ReflectionCoefficientsLocBulk RefC("../data/test_files/GreensTensorPlate.ini");
      Permittivity *perm = new PermittivityDrude("../data/test_files/GreensTensorPlate.ini");

    REQUIRE(Approx(RefC.get_epsilon(omega).real()).epsilon(1E-6) == perm->epsilon(omega).real());
  }
};
TEST_CASE("Near-field Limit of r_p works", "[ReflectionCoefficients]"){
    auto omega = GENERATE(take(5, random(-1e0, 1e0)));
    std::complex<double> kappa = 1000.;
    std::complex<double> rp, rs, rp_app, rs_app;
    ReflectionCoefficientsLocBulk RefC("../data/test_files/GreensTensorPlate.ini");
    Permittivity *perm = new PermittivityDrude("../data/test_files/GreensTensorPlate.ini");

    RefC.ref(rp , rs , omega, kappa);
      rp_app = (perm->epsilon(omega) - 1.)/(perm->epsilon(omega) + 1.);
      rs_app = (perm->epsilon(omega) - 1.)*0.25*omega*omega/(kappa*kappa);
    REQUIRE(Approx(rp.real()).epsilon(1E-4) == rp_app.real());
    REQUIRE(Approx(rp.imag()).epsilon(1E-4) == rp_app.imag());
};
TEST_CASE("r_p and r_s obey the crossing relation", "[ReflectionCoefficients]"){
    auto omega = GENERATE(take(5, random(-1e2, 1e2)));
    auto kappa_double = GENERATE(take(5, random(-1e2, 1e2)));
    std::complex<double> kappa;
    if (kappa_double < 0.){
    kappa = std::complex<double>(0.,-std::abs(kappa_double));
    };
    std::complex<double> rp_lhs, rs_lhs, rp_rhs, rs_rhs;
    ReflectionCoefficientsLocBulk RefC("../data/test_files/GreensTensorPlate.ini");
    Permittivity *perm = new PermittivityDrude("../data/test_files/GreensTensorPlate.ini");

    RefC.ref(rp_lhs , rs_lhs , omega, kappa);
    RefC.ref(rp_rhs , rs_rhs , -omega, kappa);

    REQUIRE(rp_lhs.real() == rp_rhs.real());
    REQUIRE(rs_lhs.imag() == -rs_rhs.imag());
};
