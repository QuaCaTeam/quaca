# Permittivity
```cpp
class Permittivity {
public:
  // Return the permittivity given at frequency omega
  virtual std::complex<double> epsilon(double omega) = 0;
};
```

# PermittivityDrude
```cpp
class PermittivityDrude : public Permittivity {
private:
  double omega_p, gamma; // plasma frequency and damping coefficient

public:
  // constructors
  PermittivityDrude(double omega_p, double gamma);
  PermittivityDrude(std::string input_file);

  // calculate the permittivity
  std::complex<double> epsilon(double omega);

  // getter functions
  double get_gamma() const { return this->gamma; };
  double get_omega_p() const { return this->omega_p; };
};
```
