!> TODO: Describe the tensors better. Fill in description of function. Write examples.

## GreensTensor
This is an abstract class that defines a Green's tensor.
A specific kind of Green's tensor will be a child of this class.

```cpp
class GreensTensor {
protected:
  double v;    // velocity of the particle
  double beta; // inverse temperature

public:
  // constructors
  GreensTensor(std::string input_file);
  GreensTensor(double v, double beta);

  // calculate the whole Green's tensor
  virtual void calculate_tensor(cx_mat::fixed<3, 3> &GT,
                                Options_GreensTensor opts) = 0;

  // integrate over a two-dimensional k space
  virtual void integrate_2d_k(cx_mat::fixed<3, 3> &GT,
                              Options_GreensTensor opts) = 0;

  // integrate over a one-dimensional k space
  virtual void integrate_1d_k(cx_mat::fixed<3, 3> &GT,
                              Options_GreensTensor opts) = 0;

  // getter functions
  double get_v() const { return this->v; }
  double get_beta() const { return this->beta; }

  // setter functions
  void set_v(double v) { this->v = v; }
  void set_beta(double beta) { this->beta = beta; }
};
```

### `# GreensTensor(std::string input_file)`
Input file constructor of the class.

### `# GreensTensor(double v, double beta)`
Direct constructor of the class.

### `# virtual void calculate_tensor(cx_mat::fixed<3, 3> &GT, Options_GreensTensor opts) = 0`
Calculates the Green's tensor, i.e. $\underline{G}(k, \omega + kv)$ and puts the result into the matrix `GT`.

### `# virtual void integrate_2d_k(cx_mat::fixed<3, 3> &GT, Options_GreensTensor opts) = 0`

### `# virtual void integrate_1d_k(cx_mat::fixed<3, 3> &GT, Options_GreensTensor opts) = 0`


## GreensTensorVacuum
Implements the vacuum Green's tensor given by
$$
\underline{G}_0 =
$$
where ... .

## GreensTensorPlate
Implements the Green's tensor of a plate given by
$$
\underline{G} =
$$


## Input file
The input file sections for the Green's tensor look like this
<!-- tabs:start -->
#### **GreensTensorVacuum**
```ini
[GreensTensor]
type = "vacuum"
v =
beta =
```


#### **GreensTensorPlate**
```ini
[GreensTensor]
type = "plate"
v =
beta =
```
For the plate Green's tensor you also need to define a [Permittivity](api/permittivity)!
<!-- tabs:end -->

## Examples
