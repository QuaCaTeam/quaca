//
// Created by hermasim on 03/03/2020.
//

#include "PolarizabilityNoBathMagnetic.h"

PolarizabilityNoBath::PolarizabilityNoBath(double omega_a, double alpha_zero,
                                           GreensTensor *greens_tensor)
        : PolarizabilityNoBath(omega_a, alpha_zero, greens_tensor){};

PolarizabilityNoBath::PolarizabilityNoBath(std::string input_file)
        : PolarizabilityNoBath(input_file){};

void PolarizabilityNoBath::calculate_tensor(cx_mat::fixed<3, 3> &alpha,
                                            Options_Polarizability opts) {
    // imaginary unit
    std::complex<double> I(0.0, 1.0);

    double omega = opts.omega;

    // calculate diagonal entries
    cx_mat::fixed<3, 3> diag;
    diag.zeros();
    diag(0, 0) = diag(1, 1) = diag(2, 2) = omega_a * omega_a - omega * omega;

    // calculate integral over green's tensor with fancy R
    cx_mat::fixed<3, 3> greens_R;
    struct Options_GreensTensorMagnetic opts;
    opts.omega = omega;
    opts.class_pt = this->greens_tensor;

    opts.fancy_complex = RE;
    opts.BE = RE;

    this->greens_tensor->integrate_k(greens_R, opts);

    // calculate integral over green's tensor with fancy I
    cx_mat::fixed<3, 3> greens_I;
    opts.fancy_complex = IM;
    otps.BE = IM;

    this->greens_tensor->integrate_k(greens_I, opts);

    // put everything together
    alpha =
            alpha_zero * omega_a * omega_a *
            inv(diag - alpha_zero * omega_a * omega_a * (greens_R + I * greens_I));

    if (opts.fancy_complex == IM) {
        alpha = (alpha - trans(alpha)) /
                (2.0 * I); // trans is hermitean conjugation in armadillo
    } else if (opts.fancy_complex == RE) {
        alpha = (alpha + trans(alpha)) /
                (2.0); // trans is hermitean conjugation in armadillo
    }
};
