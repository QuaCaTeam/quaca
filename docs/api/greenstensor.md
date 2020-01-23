# Green's Tensor
```cpp
class GreensTensor
{
protected:

  double v;    // velocity of the particle
  double beta; // inverse temperature

public:

  // constructors
  GreensTensor(double v, double beta);
  GreensTensor(std::string input_file);

  // calculate the whole Green's tensor
  virtual void calculate_tensor(cx_mat::fixed<3,3>& GT, vec::fixed<2> kvec, double omega) =0;

  // integrate over a two-dimensional k space
  virtual void integrate_2d_k(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts)  =0;

  // integrate over a one-dimensional k space
  virtual void integrate_1d_k(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts) =0;

  // getter functions
  double get_v();
  double get_beta();
```

# Options_GreensTensor

```cpp
struct Options_GreensTensor {
  // Different options for the integrand
  bool fancy_R = false;
  bool fancy_I = false;
  bool fancy_I_kv = false;
  bool fancy_I_temp = false;
  bool fancy_I_kv_temp = false;

  // Indices of the 3x3 GreensTensor
  vec::fixed<2> indices = {NAN, NAN};

  // Value of omega for the integration of the k-Variables
  double omega = NAN;

  // k-vector for the omega integration
  vec::fixed<2> kvec = {NAN, NAN};

  // Pointer to the GreensTensor to be able to access the attributes of the
  // class eventhough the integrand is static
  GreensTensor *class_pt;
};
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

# GreensTensorPlate
```cpp
class GreensTensorPlate : public GreensTensor
{
private:
  double z_a; // distance from plate

public:

  // constructors
  GreensTensorPlate(std::string input_file);
  GreensTensorPlate(double v, double z_a, double beta);

  // calculate full tensor
  void calculate_tensor(cx_mat::fixed<3,3>& GT, vec::fixed<2> kvec, double omega);

  // integrate over a two-dimensional k space
  void integrate_2d_k(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts);

  // integrate over a one-dimensional k space
  void integrate_1d_k(cx_mat::fixed<3,3>& GT, Options_GreensTensor opts);
};
```
