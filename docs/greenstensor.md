# Green's Tensor
```cpp
class GreensTensor
{
protected:

  double v;    // velocity of the particle
  double za;   // distance of the particle
  double beta; // inverse temperature

public:

  // constructors
  GreensTensor(double v, double za, double beta);
  GreensTensor(std::string input_file);

  // calculate the whole Green's tensor
  virtual void calculate_tensor(cx_mat::fixed<3,3>& GT, vec::fixed<2> kvec, double omega) =0;

  // integrate over a two-dimensional k space
  virtual void integrate_2d_k(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts)  =0;

  // integrate over a one-dimensional k space
  virtual void integrate_1d_k(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts) =0;

  // getter functions
  double get_v();
  double get_za();
  double get_beta();
```

# GreensTensorVacuum
```cpp
class GreensTensorVacuum : public GreensTensor
{
public:

  // constructors
  GreensTensorVacuum(double v, double beta);
  GreensTensorVacuum(std::string input_file);

  // calculate the whole Green's tensor
  void calculate_tensor(cx_mat::fixed<3,3>& GT, vec::fixed<2> kvec, double omega);

  // integrate over a two-dimensional k space
  void integrate_2d_k(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts);

  // integrate over a one-dimensional k space
  void integrate_1d_k(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts);

  // integrand for integration over one-dimensional k space
  static double integrand_1d_k(double k, void* opts);
};
```
