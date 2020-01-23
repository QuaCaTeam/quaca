# MemoryKernel

```cpp
class MemoryKernel {
public:
  // Returns the memory kernel given a frequency omega.
  virtual std::complex<double> mu(double omega) = 0;
};
```

# OhmicMemoryKernel
```cpp
class OhmicMemoryKernel : public MemoryKernel {
private:
  double gamma; // damping coefficient

public:
  // constructors
  OhmicMemoryKernel(double gamma);
  OhmicMemoryKernel(std::string input_file);

  // calculate function
  std::complex<double> mu(double omega);

  // getter functions
  double get_gamma() { return this->gamma; };
};
```
