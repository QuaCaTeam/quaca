#include "Quaca.h"
#include "catch.hpp"
#include <armadillo>
#include <complex>
#include <iomanip> // std::setprecision

TEST_CASE("The tensors from calculate_tensor and integrand_2d_k coincide",
          "[GreensTensorPlateMagnetic]") {
    // Here we considered also the volume element from the integration.
    std::complex<double> I(0.0, 1.0);
    double phi = 0;
    double k = 3.;
    double omega = 2.;
    double k_x = k*cos(phi);
    double k_y = k*sin(phi);
    double omega_p = 9;
    double gamma = 0.1;
    double v = 1e-2;
    double za = 0.1;
    double delta_cut = 30;
    vec::fixed<2> rel_err = {1E-8, 1E-6};
    double kappa_double;
    std::complex<double> kappa, volume_element;
    PermittivityDrude perm(omega_p, gamma);
    ReflectionCoefficientsLocBulk refl(&perm);
    GreensTensorPlateMagnetic Greens(v, za, 0.1, &refl, delta_cut, rel_err);
    struct Options_GreensTensorMagnetic opts;
    opts.class_pt = &Greens;

    // First, the calculate_tensor operation is used to generate the
    // Green's tensor with fancy_I
    cx_mat::fixed<3, 3> Green(fill::zeros);
    cx_mat::fixed<3, 3> Green_fancy_I(fill::zeros);
    cx_mat::fixed<3, 3> Green_fancy_R(fill::zeros);

    opts.kvec(0) = k_x;
    opts.kvec(1) = k_y;
    opts.omega = omega + k_x * v;
    k = sqrt(k_x * k_x + k_y * k_y);
    kappa = sqrt(std::complex<double>(k * k - opts.omega * opts.omega, 0.));
    kappa = std::complex<double>(std::abs(kappa.real()), -std::abs(kappa.imag()));
    volume_element = kappa * k / (k - cos(phi) * v * opts.omega);

    Greens.calculate_tensor(Green, opts);
    if (opts.omega < 0) {
        volume_element = conj(volume_element);
    }
    Green *= volume_element;
    Green_fancy_I = (Green - trans(Green)) / (2. * I);
    Green_fancy_R = (Green + trans(Green)) / (2.);

    cx_mat::fixed<3,3> Sigma = {{0,0,0},{sin(phi),-cos(phi),0},{I*kappa/k,0,-cos(phi)}};
    Sigma *= k*v/(opts.omega);

    cx_mat::fixed<3, 3> LHS(fill::zeros);
    cx_mat::fixed<3, 3> RHS(fill::zeros);

    if (kappa.real() == 0.) {
        kappa_double = kappa.imag();
    } else {
        kappa_double = kappa.real();
    }

    opts.kvec(0) = phi ;
    opts.omega = omega;
    // Second, the integrand_2d_k operation is used
    SECTION("Options: EE, IM") {
        opts.fancy_complex = IM;
        // loop over all indices
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                opts.indices(0) = i;
                opts.indices(1) = j;
                LHS(i, j) =
                        (2 * M_PI) * Greens.integrand_2d_k_magnetic(kappa_double, &opts);
                if (i != j) {
                    // As the prefactor I can not be evaluated in a purely real integration
                    // routine, it was dropped in integrate_2d_k and has to be inserted here
                    LHS(i, j) *= I;
                }
            }
        }
        RHS = Green_fancy_I;
        std::cout << LHS << std::endl << RHS << std::endl;
        std::cout << LHS/RHS << std::endl;
        REQUIRE(approx_equal(LHS, RHS, "abs", 10E-12));
    }
    SECTION("Options: EE, RE") {
        opts.fancy_complex = RE;
        // loop over all indices
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                opts.indices(0) = i;
                opts.indices(1) = j;
                LHS(i, j) =
                        (2 * M_PI) * Greens.integrand_2d_k_magnetic(kappa_double, &opts);
                if (i != j) {
                    // As the prefactor I can not be evaluated in a purely real integration
                    // routine, it was dropped in integrate_2d_k and has to be inserted here
                    LHS(i, j) *= I;
                }
            }
        }
        RHS = Green_fancy_R;
        REQUIRE(approx_equal(LHS, RHS, "abs", 10E-12));
    }
    opts.fancy_complex = IGNORE;

    SECTION("Options: BE, RE") {
        opts.BE = RE;
        // loop over all indices
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                opts.indices(0) = i;
                opts.indices(1) = j;
                LHS(i, j) =
                        (2 * M_PI) * Greens.integrand_2d_k_magnetic(kappa_double, &opts);
                if (i != j) {
                    // As the prefactor I can not be evaluated in a purely real integration
                    // routine, it was dropped in integrate_2d_k and has to be inserted here
                    LHS(i, j) *= I;
                }
            }
        }
        RHS = Sigma*Green_fancy_R;
        REQUIRE(approx_equal(LHS, RHS, "abs", 10E-12));
    }
    SECTION("Options: BE, IM") {
        opts.BE = IM;
        // loop over all indices
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                opts.indices(0) = i;
                opts.indices(1) = j;
                LHS(i, j) =
                        (2 * M_PI) * Greens.integrand_2d_k_magnetic(kappa_double, &opts);
                if (i != j) {
                    // As the prefactor I can not be evaluated in a purely real integration
                    // routine, it was dropped in integrate_2d_k and has to be inserted here
                    LHS(i, j) *= I;
                }
            }
        }
        RHS = Sigma*Green_fancy_I;
        REQUIRE(approx_equal(LHS, RHS, "abs", 10E-12));
    }
    opts.BE = IGNORE;

    SECTION("Options: EB, RE") {
        opts.EB = RE;
        // loop over all indices
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                opts.indices(0) = i;
                opts.indices(1) = j;
                LHS(i, j) =
                        (2 * M_PI) * Greens.integrand_2d_k_magnetic(kappa_double, &opts);
                if (i != j) {
                    // As the prefactor I can not be evaluated in a purely real integration
                    // routine, it was dropped in integrate_2d_k and has to be inserted here
                    LHS(i, j) *= I;
                }
            }
        }
        RHS = Green_fancy_R*trans(Sigma);
        REQUIRE(approx_equal(LHS, RHS, "abs", 10E-12));
    }
    SECTION("Options: EB, IM") {
        opts.EB = IM;
        // loop over all indices
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                opts.indices(0) = i;
                opts.indices(1) = j;
                LHS(i, j) =
                        (2 * M_PI) * Greens.integrand_2d_k_magnetic(kappa_double, &opts);
                if (i != j) {
                    // As the prefactor I can not be evaluated in a purely real integration
                    // routine, it was dropped in integrate_2d_k and has to be inserted here
                    LHS(i, j) *= I;
                }
            }
        }
        RHS = Green_fancy_I*trans(Sigma);
        REQUIRE(approx_equal(LHS, RHS, "abs", 10E-12));
    }
    opts.EB = IGNORE;

    SECTION("Options: BB, RE") {
        opts.BB = RE;
        // loop over all indices
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                opts.indices(0) = i;
                opts.indices(1) = j;
                LHS(i, j) =
                        (2 * M_PI) * Greens.integrand_2d_k_magnetic(kappa_double, &opts);
                if (i != j) {
                    // As the prefactor I can not be evaluated in a purely real integration
                    // routine, it was dropped in integrate_2d_k and has to be inserted here
                    LHS(i, j) *= I;
                }
            }
        }
        RHS = Sigma*Green_fancy_R*trans(Sigma);
        REQUIRE(approx_equal(LHS, RHS, "abs", 10E-12));
    }
    SECTION("Options: BB, IM") {
        opts.BB = IM;
        // loop over all indices
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                opts.indices(0) = i;
                opts.indices(1) = j;
                LHS(i, j) =
                        (2 * M_PI) * Greens.integrand_2d_k_magnetic(kappa_double, &opts);
                if (i != j) {
                    // As the prefactor I can not be evaluated in a purely real integration
                    // routine, it was dropped in integrate_2d_k and has to be inserted here
                    LHS(i, j) *= I;
                }
            }
        }
        RHS = Sigma*Green_fancy_I*trans(Sigma);
        REQUIRE(approx_equal(LHS, RHS, "abs", 10E-12));
    }
    opts.BB = IGNORE;
};

