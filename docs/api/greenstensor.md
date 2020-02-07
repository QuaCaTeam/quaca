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
Implements the imaginary part of the vacuum Green's tensor given by
$$
  \mathrm{Im}\left\{\underline{G}_0(\mathbf{ k}, z- z'\to 0,\omega)\right\} =
  \frac{
  \theta(\frac{\omega^2}{c^2}-k^2)
}{2\epsilon_0\sqrt{\omega^2/c^2-k^2}}
  \left(
    \mathbb{1}\frac{\omega^2}{c^2} -
    \mathrm{diag}\left[
      k_x^2,\,k_y^2,\,\frac{\omega^2}{c^2}-k^2
    \right]
  \right)
  ,
$$
where $\omega$ is the frequency, $\mathbf{k}=(k_x,\,k_y)^\intercal$ is the corresponding wavevector to the two-dimensional $xy$ plane (with $|\mathbf{k}|=k$, and $z$ is the remaining (not Fourier-transformed) spatial coordinate. The $z- z'\to 0$ indicates, that the Green's tensor is evaluated twice at the same point with respect to the $z$ coordinate. The introduced Heaviside function $\theta(\omega^2/c^2-k^2)$ ensures, that solely imaginary values are considered. The real part of the vacuum Green's tensor is implicitly contained in the $\omega_a$. A clear derivation of the vacuum Green's tensor can be found for instance in [this paper](https://www.mdpi.com/2076-3417/7/11/1158).
## GreensTensorPlate
Implements the scattered part of the Green's tensor of a flat surface beneath vaccuum and is given by
$$
  \underline{g}(\mathbf{k},z_a, z_a,\omega) =
  \frac{\omega^2}{c^2}
  \big(
  r^s(\mathbf{k},\omega) \, \mathbf{ e}_s \otimes {\mathbf{ e}_s} 
  +
  r^p(\mathbf{k},\omega) \,\cdot \mathbf{ e}_{+p} \otimes {\mathbf{ e}_{-p}} 
\big)\frac{e^{-2\kappa z_a}}{2\epsilon_0\kappa}
,
$$
where $r^s$ and $r^p$ are the reflection coefficients of the $s$ (transverse electric) and $p$ (transverse magnetic) polarization. The corresponding unit vector are
$$
  \mathbf{e}_s = \frac{\mathbf{k}}{k}\times \frac{\mathbf{z}}{z}  = \big(
    -\frac{k_y}{k},\, \frac{k_x}{k},\, 0
\big)^\intercal
\quad \text{and} \quad
\mathbf{ e}_{\pm p} = \frac{c\kappa}{\omega} \left(\frac{\mathbf{ z}}{z}\frac{k}{\kappa} \mp \mathrm{i} \frac{\mathbf{ k}}{k}\right) =\frac{c \kappa}{\omega} \big(
  \mp\mathrm{i}\frac{k_x}{k},\,\mp\mathrm{i}\frac{k_y}{k},\, \frac{k}{\kappa}.
\big)^\intercal.
$$
A reference for the scattered part of the Green's tensor can be found in [this paper](http://link.aps.org/doi/10.1103/PhysRevA.94.042114). Notice that the full Green's tensor of the described configuration reads $\underline{G}=\underline{G}_0+\underline{g}$. 
## Input file
The input file sections for the Green's tensor look like this
<!-- tabs:start -->
#### **GreensTensorVacuum**
```ini
[GreensTensor]
type = vacuum
v =
beta =
```


#### **GreensTensorPlate**
```ini
[GreensTensor]
type = plate
v =
beta =
```
For the plate Green's tensor you also need to define the [Reflection Coefficients](api/reflection)!
<!-- tabs:end -->

## Examples
