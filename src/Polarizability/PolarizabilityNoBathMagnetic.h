//
// Created by hermasim on 03/03/2020.
//

#ifndef QUACA_POLARIZABILITYNOBATHMAGNETIC_H
#define QUACA_POLARIZABILITYNOBATHMAGNETIC_H

#include "../GreensTensor/GreensTensor.h"
#include "../GreensTensor/GreensTensorPlateMagnetic.h"
#include "Polarizability.h"
#include "PolarizabilityNoBath.h"


class PolarizabilityNoBathMagnetic: public PolarizabilityNoBath {
public:
  GreensTensorPlateMagnetic* greens_magnetic;
    PolarizabilityNoBathMagnetic(double omega_a, double alpha_zero,
                         GreensTensor *greens_tensor);
    PolarizabilityNoBathMagnetic(std::string input_file);

    void calculate_tensor(cx_mat::fixed<3, 3> &alpha,
                          Options_Polarizability opts);
};


#endif //QUACA_POLARIZABILITYNOBATHMAGNETIC_H
