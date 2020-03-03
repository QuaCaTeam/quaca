//
// Created by hermasim on 20/02/2020.
//
#include "GreensTensorPlateMagnetic.h"
#include <cmath>
#include <complex>


void GreensTensorPlateMagnetic::calculate_tensor(cx_mat::fixed<3, 3> &GT,
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
    GT(1,0) = (-prefactor_s*ky*kx/k_quad + prefactor_p*ky*kx/k_quad)*kappa;
    GT(0,1) = GT(1,0);
    GT(2,1) = prefactor_p*I*ky;
    GT(1,2) = - GT(2,1);

    // In case of negative frequencies, the tensor has to be hermitian transposed
    if (opts.omega < 0) {
        GT = trans(GT);
    }
};

double GreensTensorPlateMagnetic::integrand_2d_k_magnetic(double kappa_double, void *opts) {
    // The needed parameters for the integration are encoded in the void pointer.
    // This void pointer is casted the options struct given in GreensTensor.h.
    Options_GreensTensorMagnetic *opts_pt = static_cast<Options_GreensTensorMagnetic *>(opts);
    // To access attributes of the GreensTensorPlate class, the class pointer
    // within the struct is casted.
    GreensTensorPlateMagnetic *pt = static_cast<GreensTensorPlateMagnetic *>(opts_pt->class_pt);

    // read general input parameters
    double beta = pt->beta;
    double v = pt->v;
    double v_quad = v * v;
    double za = pt->za;
    // read omega and phi from the option struct
    double omega = opts_pt->omega;
    double omega_quad = omega * omega;
    double phi = opts_pt ->kvec(0);
    double cos_quad = pow(cos(phi),2);
    double sin_quad = 1.0 - cos_quad;

    // Before the real or imaginary part of the chosen matrix element can be
    // calculated, the complex result is stored in result_complex.
    std::complex<double> result_complex = (0,0);
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

    // Transfer kappa to the correct complex value
    if (kappa_double < 0.0) {
        kappa_complex = std::complex<double>(0.0, kappa_double);
        kappa_quad = -kappa_double * kappa_double;
    } else {
        kappa_complex = std::complex<double>(kappa_double, 0.0);
        kappa_quad = kappa_double * kappa_double;
    }

    // Express k via frequency and kappa
    k = (sqrt(kappa_quad * (1.0 - v_quad * cos_quad) + omega_quad) +
         v * omega * cos(phi)) /
        (1 - v_quad * cos_quad);
    k_quad = k * k;

    // Define the Doppler-shifted frequency
    omega_pl = (omega + k * cos(phi) * v);
    omega_pl_quad = omega_pl * omega_pl;

    // In order to obey reality in time, a positive omega_pl is used for the
    // actual calculation. Afterwards, the corresponding symmetry operation is
    // performed if the sign of omega_pl is negative.
    double omega_pl_abs = std::abs(omega_pl);

    // producing the reflection coefficients in p- and s-polarization
    pt->reflection_coefficients->ref(r_p, r_s, omega_pl_abs, kappa_complex);

    //Enforce the crossing relation
    if(omega_pl < 0)
    {
        kappa_complex = conj(kappa_complex);
        r_s = conj(r_s);
        r_p = conj(r_p);
    }

    //check wether the electric Green's tensor should be computed
    if(opts_pt->fancy_complex != IGNORE) {

        //General prefactor for all components
        std::complex<double> prefactor = kappa_complex/(1.-cos(phi)*v*omega_pl/k);
        //prefactors for the ss-polarization, the hermitian part of the pp polarization
        //and the rotational part of the pp-polarization
        std::complex<double> r_s_factor;
        std::complex<double> r_p_herm_factor;
        std::complex<double> r_p_rot_factor;
        if(opts_pt->fancy_complex == IM)
        {
            r_s_factor = prefactor*omega_pl_quad*imag(r_s/kappa_complex*exp(-2*za*kappa_complex));
            r_p_herm_factor = prefactor*imag(r_p*kappa_complex*exp(-2.*za*kappa_complex));
            r_p_rot_factor = prefactor*kappa_complex*imag(r_p*exp(-2.*za*kappa_complex));
        }
        if(opts_pt->fancy_complex == RE)
        {
            r_s_factor = prefactor*omega_pl_quad*real(r_s/kappa_complex*exp(-2*za*kappa_complex));
            r_p_herm_factor = prefactor*real(r_p*kappa_complex*exp(-2.*za*kappa_complex));
            r_p_rot_factor = prefactor*kappa_complex*real(r_p*exp(-2.*za*kappa_complex));
        }

        // Calculate the G_xx element
        if (opts_pt->indices(0) == 0 && opts_pt->indices(1) == 0) {
            result_complex +=  r_s_factor * sin_quad;
            result_complex += r_p_herm_factor * cos_quad;
        }
            // Calculate the G_yy element
        else if (opts_pt->indices(0) == 1 && opts_pt->indices(1) == 1) {
            result_complex += r_p_herm_factor * sin_quad;
            result_complex += r_s_factor * cos_quad;
        }
            // Calculate the G_zz element
        else if (opts_pt->indices(0) == 2 && opts_pt->indices(1) == 2) {
            result_complex += r_p_herm_factor * k_quad / kappa_quad;
        }
            // Calculate the G_zx element
        else if (opts_pt->indices(0) == 2 && opts_pt->indices(1) == 0) {
            result_complex += r_p_rot_factor*I*k*cos(phi)/kappa_complex;
        }
            // Calculate the G_xz element
        else if (opts_pt->indices(0) == 0 && opts_pt->indices(1) == 2) {
            result_complex += -r_p_rot_factor* I * k*cos(phi)/kappa_complex;
        }
    }
    //reset the result_complex variable to store additonal components of the other Green's tensors
    //check if the BE-GreensTensor should be computed
    if(opts_pt->BE != IGNORE)
    {
        //General prefactor for all components
        std::complex<double> prefactor = kappa_complex/(1.-cos(phi)*v*omega_pl/k);
        //prefactors for the ss-polarization, the hermitian part of the pp polarization
        //and the rotational part of the pp-polarization
        std::complex<double> r_s_factor;
        std::complex<double> r_p_herm_factor;
        std::complex<double> r_p_rot_factor;
        if(opts_pt-> BE == IM)
        {
            r_s_factor = prefactor*omega_pl_quad*imag(r_s/kappa_complex*exp(-2*za*kappa_complex));
            r_p_herm_factor = prefactor*imag(r_p*kappa_complex*exp(-2.*za*kappa_complex));
            r_p_rot_factor = prefactor*kappa_complex*imag(r_p*exp(-2.*za*kappa_complex));
        }
        if(opts_pt-> BE == RE)
        {
            r_s_factor = prefactor*omega_pl_quad*real(r_s/kappa_complex*exp(-2*za*kappa_complex));
            r_p_herm_factor = prefactor*real(r_p*kappa_complex*exp(-2.*za*kappa_complex));
            r_p_rot_factor = prefactor*kappa_complex*real(r_p*exp(-2.*za*kappa_complex));
        }
       //Compute all the different components of the Green's tensor
       //G_zx
       if(opts_pt->indices(0) == 2 && opts_pt->indices(1) == 0)
       {
           result_complex += r_s_factor*I*v*kappa_complex*sin_quad/omega_pl;
           result_complex += r_p_herm_factor*I*v*kappa_complex*cos_quad/omega_pl;
           result_complex -=r_p_rot_factor*I*k_quad*v*cos_quad/(omega_pl*kappa_complex);
       }
       //G_yy
       else if(opts_pt->indices(0) == 1 && opts_pt->indices(1) == 1)
       {
           result_complex -= r_s_factor*k*v*cos(phi)/omega_pl;
       }
       //G_zz
       else if(opts_pt->indices(0) == 2 && opts_pt->indices(1) == 2)
       {
           result_complex -= r_p_herm_factor*pow(k,3)*v*cos(phi)/(kappa_quad*omega_pl);
           result_complex += r_p_rot_factor*k*v*cos(phi)/(omega_pl);
       }
    }
    if(opts_pt->EB != IGNORE)
    {
        //General prefactor for all components
        std::complex<double> prefactor = kappa_complex/(1-cos(phi)*v*omega_pl/k);
        //prefactors for the ss-polarization, the hermitian part of the pp polarization
        //and the rotational part of the pp-polarization
        std::complex<double> r_s_factor;
        std::complex<double> r_p_herm_factor;
        std::complex<double> r_p_rot_factor;
        if(opts_pt-> EB == IM)
        {
            r_s_factor = prefactor*omega_pl_quad*imag(r_s/kappa_complex*exp(-2*za*kappa_complex));
            r_p_herm_factor = prefactor*imag(r_p*kappa_complex*exp(-2.*za*kappa_complex));
            r_p_rot_factor = prefactor*kappa_complex*imag(r_p*exp(-2.*za*kappa_complex));
        }
        if(opts_pt-> EB == RE)
        {
            r_s_factor = prefactor*omega_pl_quad*real(r_s/kappa_complex*exp(-2*za*kappa_complex));
            r_p_herm_factor = prefactor*real(r_p*kappa_complex*exp(-2.*za*kappa_complex));
            r_p_rot_factor = prefactor*kappa_complex*real(r_p*exp(-2.*za*kappa_complex));
        }
        //Compute all the different components of the Green's tensor
        //G_xz
        if(opts_pt->indices(0) == 0 && opts_pt->indices(1) == 2)
        {
            result_complex -= r_s_factor*I*v*conj(kappa_complex)*sin_quad/omega_pl;
            result_complex -= r_p_herm_factor*I*v*conj(kappa_complex)*cos_quad/omega_pl;
            result_complex +=r_p_rot_factor*I*k_quad*v*cos_quad/(omega_pl*kappa_complex);
        }
        //G_yy
        else if(opts_pt->indices(0) == 1 && opts_pt->indices(1) == 1)
        {
            result_complex -= r_s_factor*k*v*cos(phi)/omega_pl;
        }
        //G_zz
        else if(opts_pt->indices(0) == 2 && opts_pt->indices(1) == 2)
        {
            result_complex -= r_p_herm_factor*pow(k,3)*v*cos(phi)/(kappa_quad*omega_pl);
            result_complex += r_p_rot_factor*k*v*cos(phi)*conj(kappa_complex)/(kappa_complex*omega_pl);
        }

    }
    if(opts_pt->BB != IGNORE)
    {
        //General prefactor for all components
        std::complex<double> prefactor = kappa_complex/(1-cos(phi)*v*omega_pl/k);
        //prefactors for the ss-polarization, the hermitian part of the pp polarization
        //and the rotational part of the pp-polarization
        std::complex<double> r_s_factor;
        std::complex<double> r_p_herm_factor;
        std::complex<double> r_p_rot_factor;
        if(opts_pt-> BB == IM)
        {
            r_s_factor = prefactor*omega_pl_quad*imag(r_s/kappa_complex*exp(-2*za*kappa_complex));
            r_p_herm_factor = prefactor*imag(r_p*kappa_complex*exp(-2.*za*kappa_complex));
            r_p_rot_factor = prefactor*kappa_complex*imag(r_p*exp(-2.*za*kappa_complex));
        }
        if(opts_pt-> BB == RE)
        {
            r_s_factor = prefactor*omega_pl_quad*real(r_s/kappa_complex*exp(-2*za*kappa_complex));
            r_p_herm_factor = prefactor*real(r_p*kappa_complex*exp(-2.*za*kappa_complex));
            r_p_rot_factor = prefactor*kappa_complex*real(r_p*exp(-2.*za*kappa_complex));
        }
        //Compute all the different components of the Green's tensor
        //G_yy
        if(opts_pt->indices(0) == 1 && opts_pt->indices(1) == 1)
        {
            result_complex += r_s_factor*k_quad*pow(v,2)/omega_pl_quad;
        }
        //G_zz
        else if(opts_pt->indices(0) == 2 && opts_pt->indices(1) == 2)
        {
            result_complex += r_s_factor*pow(v,2)*kappa_quad*sin_quad/omega_pl_quad;
            result_complex += r_p_herm_factor*pow(v,2)*(pow(k,4)+pow(kappa_complex,3)*conj(kappa_complex))
                    *cos_quad/(kappa_quad*omega_pl_quad);
            result_complex -= r_p_rot_factor*2.*k_quad*pow(v,2)*cos_quad*real(kappa_complex)
                    /(kappa_complex*omega_pl_quad);
        }
    }
    //Depending on the component that was computed either the real or the imaginary part should be zero. Therefore
    //we can add both of the up. The missing I has to be added after integration procedure
    if(real(result_complex) == 0)
    {
        result = imag(result_complex);
    }
    else if(imag(result_complex) == 0)
    {
        result = real(result_complex);
    }
    else {
        std::cerr << "Result_complex is neither purely real nor purely complex. Such a result is not implemented yet" << std::endl;
        exit(0);
    }

    // Add weighting function if demanded
    if (opts_pt->weight_function == KV) {
        result *= k * cos(phi);
    } else if (opts_pt->weight_function == TEMP) {
        result /= (1.0 - exp(-beta * omega_pl));
    } else if (opts_pt->weight_function == NON_LTE) {
        result *=
                1. / (1.0 - exp(-beta * omega_pl)) - 1. / (1.0 - exp(-beta * omega));
    } else if (opts_pt->weight_function == KV_TEMP) {
        result *= k * cos(phi) / (1.0 - exp(-beta * omega_pl));
    } else if (opts_pt->weight_function == KV_NON_LTE) {
        result *=
                k * cos(phi) *
                (1. / (1.0 - exp(-beta * omega_pl)) - 1. / (1.0 - exp(-beta * omega)));
    }


    return result;
};

