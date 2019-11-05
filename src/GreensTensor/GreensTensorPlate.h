#ifndef GREENSTENSORPLATE_H
#define GREENSTENSORPLATE_H

class GreensTensorPlate : public GreensTensor
{
protected:
  cx_mat::fixed<3,3> greens;
  double v, za;
  Permittivity *permittivity;

public:

  GreensTensorPlate(std::string input_file);
  cx_mat::fixed<3,3> calculate_pure(double k_v, double omega);
  cx_mat::fixed<3,3> calculate_integrated(int flags, double omega);

};


#endif // GREENSTENSORPLATE_H
