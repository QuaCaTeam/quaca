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
Calculates the Green's tensor integrated along a one-dimensional k-vector, such that the result depends on the frequenccy $\omega$ and a one-dimensional k-vector $k$ e.g
$$
\int d k_1 \underline{G}(\mathbf{k},z_a,\omega) = \underline{G}(k_2,z_a,\omega)
$$
The chosen basis for the two-dimensional k-vector (e.g $(k_x,k_y)$ or $(k,\phi)$) depends on the implementation of the Green's tensor. The result is stored in the matrix `GT`.

### `# virtual void integrate_1d_k(cx_mat::fixed<3, 3> &GT, Options_GreensTensor opts) = 0`
Calculates the Green's tensor integrated along the two-dimensional k-vector, such that the result only depends on the frequency $\omega$ and the distance $z_a$. The result is stored in  the matrix `GT`.

## GreensTensorVacuum
Implements the vacuum Green's tensor given by
$$
  \underline{G}_0(\mathbf{k}, z,z',\omega) =
  \frac{\omega^2}{\epsilon_0c^2}
  \mathcal{P}\left[
    \mathbb{1}-\frac{c^2}{\omega^2}
    \begin{bmatrix}
      \mathbf{k} \\ \pm i \kappa
    \end{bmatrix}
    \left[ \mathbf{k}^\intercal ,\,\pm i\kappa\right]
  \right]
  \frac{e^{-\kappa|z-z'|}}{2\kappa}
  -\frac{\mathbf{e}_z\otimes\mathbf{e}_z}{\omega^2/c^2}\delta(z-z')
  .
$$
Here the case $+$ relates to $z > z'$ and $-$ to $z < z '$.
$\kappa=\sqrt{k^2-\omega^2/c^2}$ with $\mathrm{Re}\{\kappa\}\geq0$ and $\mathrm{Im}\{\kappa\}<0$. Especially, we focus on the case $z=z'$. As the real part of this Green's tensor is divergent for this case, we can continue calculating the imaginary part. Here we also eliminated odd orders of $k_y$, since we solely consider symmetric integration over $k_y$
$$
  \mathrm{Im}\left\{\underline{G}_0(\mathbf{k}, z\to z',\omega)\right\} =
  \frac{
  \theta(\frac{\omega^2}{c^2}-k^2)
}{2\epsilon_0\sqrt{\omega^2/c^2-k^2}}
  \left(
    \mathbb{1}\frac{\omega^2}{c^2} -
    \mathrm{diag}\left[
      k_x^2,\,k_y^2,\,\frac{\omega^2}{c^2}-k^2
    \right]
  \right)
$$ 
```cpp
class GreensTensorVacuum : public GreensTensor {
public:
  // constructors
  GreensTensorVacuum(double v, double beta);
  GreensTensorVacuum(std::string input_file);

  // calculate the tensor in frequency and momentum space
  void calculate_tensor(cx_mat::fixed<3, 3> &GT, Options_GreensTensor opts);

  // integrate over a two-dimensional k space
  void integrate_2d_k(cx_mat::fixed<3, 3> &GT, Options_GreensTensor opts);

  // integrate over a one-dimensional k space
  void integrate_1d_k(cx_mat::fixed<3, 3> &GT, Options_GreensTensor opts);

  // integrand for integration over one-dimensional k space
  static double integrand_1d_k(double k, void *opts);
};
```

###  `void calculate_tensor(cx_mat::fixed<3, 3> &GT, Options_GreensTensor opts)`
Calculates the vacuum Green's tensor as given above. As the real part is divergent only the imaginary part can be calculated. Trying to compute the real part will result in an abort of the computation and an error message.


###  `void integrate_2d_k(cx_mat::fixed<3, 3> &GT, Options_GreensTensor opts)`
Computes the vaccum Green's tensor integrated along the $k_y$ direction. The result is analytically obtained and takes the form
$$
\int d k_y \text{Im}\{\underline{G}(\mathbf{k},z\to z', \omega + \mathbf{k}\cdot\mathbf{v}) \} = \\ \text{diag}\big[ \xi^2,(\omega+k_x v)^2/c^2-\xi^2/2,(\omega+k_x v)^2/c^2-\xi^2/2\big]
$$
!> There are some prefactors missing

###  `void integrate_1d_k(cx_mat::fixed<3, 3> &GT, Options_GreensTensor opts)`
Computes the integration along the 2-dimensional $\mathbf{k}$-vector. There are different additional factors in the integrand, which can be set by the `Options_GreensTensor` struct. For more information see the description of [Options_GreensTensor](#Options_GreensTensor).

### `static double integrand_1d_k(double k, void *opts)`

## GreensTensorPlate
Implements the Green's tensor of a plate given by
$$
\underline{G} =
$$

# Options_GreensTensor

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
For the plate Green's tensor you also need to define a [Permittivity](api/permittivity)!
<!-- tabs:end -->

## Examples
