#ifndef GREENSTENSOR_H
#define GREENSTENSOR_H

// A Greens tensor class
/*!
 * This is an abstract class that implements an isotropic and reciprocal Greens tensor.
 */
class GreensTensor
{
private:
    std::complex<double> greens[3][3];
    double v, za;
};



#endif // GREENSTENSOR_H
