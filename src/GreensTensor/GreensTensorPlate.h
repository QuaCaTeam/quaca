#ifndef GREENSTENSORPLATE_H
#define GREENSTENSORPLATE_H

class GreensTensorPlate : public GreensTensor
{
protected:
  double v, za;

public:

  GreensTensorPlate(std::string input_file);
  void calculate_pure(cx_mat::fixed<3,3>& GT, vec::fixed<2> kvec, double omega);
  void calculate_integrated(cx_mat::fixed<3,3>& GT, double omega, Options *opts);

};


#endif // GREENSTENSORPLATE_H
