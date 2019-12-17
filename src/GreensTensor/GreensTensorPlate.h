#ifndef GREENSTENSORPLATE_H
#define GREENSTENSORPLATE_H

class GreensTensorPlate : public GreensTensor
{
public:

  GreensTensorPlate(std::string input_file);
  void calculate_tensor(cx_mat::fixed<3,3>& GT, vec::fixed<2> kvec, double omega);
  void integrate_k_2d(cx_mat::fixed<3,3>& GT, Options_GreensTensor *opts);
  void integrate_k_1d(cx_mat::fixed<3,3>& GT, Options_GreensTensor *opts);

};


#endif // GREENSTENSORPLATE_H
