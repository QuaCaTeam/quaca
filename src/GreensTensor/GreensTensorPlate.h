#ifndef GREENSTENSORPLATE_H
#define GREENSTENSORPLATE_H

class GreensTensorPlate : public GreensTensor
{
protected:
  cx_mat::fixed<3,3> greens;
  double v, za;

public:

  GreensTensorPlate(std::string input_file);
  void calculate_pure(cx_mat::fixed<3,3> GT, vec::fixed<2> kvec, double omega);
  void calculate_integrated(cx_mat::fixed<3,3> GT, double omega, std::vector<std::string> options);

};


#endif // GREENSTENSORPLATE_H
