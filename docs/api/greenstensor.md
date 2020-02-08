!> TODO: Write examples.

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
Calculates the Green's tensor, i.e. $\underline{G}(\mathbf{k}, \omega + k_x v)$, where $x$ is the axis along the motion with velocity $v$, and puts the result into the matrix `GT`.

### `# virtual void integrate_1d_k(cx_mat::fixed<3, 3> &GT, Options_GreensTensor opts) = 0`
Calculates the integral over the Green's tensor and writes it into the matrix `GT`, i.e. $\int \frac{\mathrm{d}^2\mathbf{k}}{(2\pi)} \underline{G}(\mathbf{k},\omega+k_xv) $. Here following integration options are available and governed by the flags in `opts`.

| Flag                | Math       |
|---------------------|------------|
| `fancy_R`           | $\underline{G}_\Re = (\underline{G}+\underline{G}^\dagger)/2$ |
| `fancy_I`           | $\underline{G}_\Im = (\underline{G}-\underline{G}^\dagger)/(2\mathrm{i})$ |
| `fancy_I_kv`           | $k_x \, \underline{G}_\Im $ |
| `fancy_I_temp`           | $ \frac{\underline{G}_\Im}{1-\exp(-\hbar\beta(\omega+k_xv))} $ |
| `fancy_I_non_LTE`           | $ \underline{G}_\Im \left(\frac{1}{1-\exp(-\hbar\beta(\omega+k_xv))} - \frac{1}{1-\exp(-\hbar\beta\omega)}\right) $ |
| `fancy_I_temp`           | $ \frac{k_x\,\underline{G}_\Im}{1-\exp(-\hbar\beta(\omega+k_xv))} $ |
| `fancy_I_non_LTE`           | $ k_x\,\underline{G}_\Im \left(\frac{1}{1-\exp(-\hbar\beta(\omega+k_xv))} - \frac{1}{1-\exp(-\hbar\beta\omega)}\right) $ |

If a flag is `true` the corresponding integrand of the above table is calculated. Setting multiple flags to `true` might not lead to the desired outcome.

### `# virtual void integrate_2d_k(cx_mat::fixed<3, 3> &GT, Options_GreensTensor opts) = 0`
Performs the first of the two integrations and analogously to `integrate_2d_k`. However, which specific variable is integrated first and second depends on the child of the class. 

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
A reference for the scattered part of the Green's tensor can be found in [this paper](http://link.aps.org/doi/10.1103/PhysRevA.94.042114). Please note, that the odd orders of $k_y$ are, due to symmetry reasons, not considered. Therefore, the elements $G_{xy}=G_{yx}=G_{yz}=G_{zy}=0$. Moreover, notice that the full Green's tensor of the described configuration reads $\underline{G}=\underline{G}_0+\underline{g}$. 
## Input file
The input file sections for the Green's tensor look like this
<!-- tabs:start -->
#### **GreensTensorVacuum**
```ini
[GreensTensor]
type = vacuum
v =
beta =
rel_err_1 =
```


#### **GreensTensorPlate**
```ini
[GreensTensor]
type = plate
v =
beta =
za = 
delta_cut =
rel_err_0 =
rel_err_1 =
```
For the plate Green's tensor you also need to define the [Reflection Coefficients](api/reflection)!
<!-- tabs:end -->

## Examples
<!-- tabs:start -->

#### ** Example: Vacuum **
We want to calculate the twofold integrated imaginary part of the vaccuum Green's tensor at frequency $\omega=3\,\mathrm{eV}$ in the static case $(v=0)$, and at $\beta=0.1\,\mathrm{eV}^{-1}$. Moreover, we need to define the accuracy of the numerical integration. Since the first integration of the vacuum Green's tensor is implemented analytically, we solely need to provide the demanded relative accuracy for the second integration $\delta_\mathrm{err}=10^{-9}$.
Let's define the Green's tensor
```cpp
double v = 0;
double beta = 1e-1;
double relerr = 1e-9
GreensTensorVacuum greens_tensor(v, beta, relerr);
```
To calculate the integration we perform we pass the desired $\omega$ and the integration option `fancy_I` to `opts`
```cpp
cx_mat::fixed<3,3> GT(fill::zeros);
Options_GreensTensor opts;
opts.omega = 3.0;
opts.fancy_I = true;
greens_tensor.integrated_1d_k(GT, opts);
```
The matrix `GT` now contains the integrated imaginary part of the vacuum Green's tensor.


#### ** Example (.ini): Vacuum **
We want to calculate the twofold integrated imaginary part of the vaccuum Green's tensor at frequency $\omega=3\,\mathrm{eV}$ in the static case $(v=0)$, and at $\beta=0.1\,\mathrm{eV}^{-1}$. Moreover, we need to define the accuracy of the numerical integration. Since the first integration of the vacuum Green's tensor is implemented analytically, we solely need to provide the demanded relative accuracy for the second integration $\delta_\mathrm{err}=10^{-9}$.
Since we have to define a lot of parameters, QuaCa offers a shortcut to the long task before.
We simply define all parameters in a file called `parameters.ini` which looks like this
```ini
[GreensTensor]
type = vacuum
v = 0
beta = 1e-1
rel_err_1 = 1e-9
```
Now we can easily define the vacuum Green's tensor class and calculate the desired part of the integrated Green's tensor
```cpp
GreensTensorVacuum greens_tensor("parameters.ini");

cx_mat::fixed<3,3> GT(fill::zeros);
Options_GreensTensor opts;
opts.omega = 3.0;
opts.fancy_I = true;
greens_tensor.integrated_1d_k(GT, opts);
```



<!-- tabs:end -->
