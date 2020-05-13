# GreensTensor {docsify-ignore-all}
This is an abstract class that defines a Green's tensor.
A specific kind of Green's tensor will be a child of this class.

```cpp
class GreensTensor {
protected:
  double v; //velocity of the particle 
  double beta; // inverse temperature

public:
  // constructor
  GreensTensor(double v, double beta);
  explicit GreensTensor(const std::string &input_file);

  // calculate the tensor in frequency and momentum space
  virtual void calculate_tensor(double omega, vec::fixed<2> k,
                                cx_mat::fixed<3, 3> &GT) const = 0;

  // integrate over a two-dimensional k space
  virtual void integrate_k(double omega, cx_mat::fixed<3, 3> &GT,
                           Tensor_Options fancy_complex,
                           Weight_Options weight_function) const = 0;

  // getter functions
  double get_v() const { return this->v; };
  double get_beta() const { return this->beta; };

  virtual double omega_ch() const = 0;

  // setter function
  virtual void set_v(double v_new) { this->v = v_new; };
};
```

### `# GreensTensor(std::string input_file)`
Input file constructor of the class.

### `# GreensTensor(double v, double beta)`
Direct constructor of the class.

### `# virtual void calculate_tensor(double omega, vec::fixed<2> k,cx_mat::fixed<3, 3> &GT) const = 0;`
Calculates the Green's tensor, i.e. $\underline{G}(\mathbf{k}, \omega + k_x v)$, where $x$ is the axis along the motion with velocity $v$, and puts the result into the matrix `GT`.

###  `# virtual void integrate_k(double omega, cx_mat::fixed<3, 3> &GT, Tensor_Options fancy_complex, Weight_Options weight_function) const = 0;`
Compute the integral over the momentum space of the Green's tensor i.e. $\int_{-\infty}^{\infty} \frac{\mathrm{d} \mathbf{k}}{(2 \pi)^2} \; f(\mathbf{k},\omega) \underline{G}(\mathbf{k}, \omega + k_x v)$, where $f(\mathbf{k},\omega)$ is some weight function.  The result is stored in `GT`.
The different options for the two variables `fancy_complex` and `weight_function` are

| Flag                | Math       |
|---------------------|------------|
| `fancy_R`           | $\underline{G}_\Re = (\underline{G}+\underline{G}^\dagger)/2$ |
| `fancy_I`           | $\underline{G}_\Im = (\underline{G}-\underline{G}^\dagger)/(2\mathrm{i})$ |
| `fancy_I_kv`           | $k_x \, \underline{G}_\Im $ |
| `fancy_I_temp`           | $ \frac{\underline{G}_\Im}{1-\exp(-\hbar\beta(\omega+k_xv))} $ |
| `fancy_I_non_LTE`           | $ \underline{G}_\Im \left(\frac{1}{1-\exp(-\hbar\beta(\omega+k_xv))} - \frac{1}{1-\exp(-\hbar\beta\omega)}\right) $ |
| `fancy_I_temp`           | $ \frac{k_x\,\underline{G}_\Im}{1-\exp(-\hbar\beta(\omega+k_xv))} $ |
| `fancy_I_non_LTE`           | $ k_x\,\underline{G}_\Im \left(\frac{1}{1-\exp(-\hbar\beta(\omega+k_xv))} - \frac{1}{1-\exp(-\hbar\beta\omega)}\right) $ |

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
```cpp

class GreensTensorVacuum : public GreensTensor {
private:
  double relerr; //integration error along the k_v direction

public:
  // constructors
  GreensTensorVacuum(double v, double beta, double relerr);
  GreensTensorVacuum(std::string input_file);

  // calculate the tensor in frequency and momentum space
  void calculate_tensor(double omega, vec::fixed<2> k,
                        cx_mat::fixed<3, 3> &GT) const;

  // integrate over a two-dimensional k space
  void integrate_k(double omega, cx_mat::fixed<3, 3> &GT,
                   Tensor_Options fancy_complex,
                   Weight_Options weight_function) const;

  // integrand for integration over one-dimensional k space
  double integrand_k(double kv, double omega, const vec::fixed<2> &indices,
                                         Tensor_Options fancy_complex,
                                         Weight_Options weight_function) const;

  double omega_ch() const;
  double get_relerr() const { return this->relerr; };
};

```
### `# double integrand_k(double kv, double omega, const vec::fixed<2> &indices,Tensor_Options fancy_complex,Weight_Options weight_function) const;`
                                         
Implements the integrand in momentum space along the direction of motion of the microscopic particle.
                                         


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

```cpp
class GreensTensorPlate : public GreensTensor {
protected:
  //distance between microscopic object and macroscopic surface
  double za;

  // kappa_cut defines the numerical cut-off of the kappa integration
  double delta_cut;
  vec::fixed<2> rel_err = {NAN, NAN};

  // reflection coefficients are needed to describe the surface's response
  std::shared_ptr<ReflectionCoefficients> reflection_coefficients;

public:
  // constructors
  GreensTensorPlate(
      double v, double beta, double za,
      std::shared_ptr<ReflectionCoefficients> reflection_coefficients,
      double delta_cut, const vec::fixed<2> &rel_err);
  explicit GreensTensorPlate(const std::string &input_file);

  // calculate the tensor in frequency and momentum space
  void calculate_tensor(double omega, vec::fixed<2> k,
                        cx_mat::fixed<3, 3> &GT) const;

  // integrate over a two-dimensional k space
  void integrate_k(double omega, cx_mat::fixed<3, 3> &GT,
                   Tensor_Options fancy_complex,
                   Weight_Options weight_function) const;

  // integrands
  double integrand_1d_k(double phi, double omega, const vec::fixed<2> &indices,
                        Tensor_Options fancy_complex,
                        Weight_Options weight_function) const;

  double integrand_2d_k(double kappa_double, double omega, double phi,
                        const vec::fixed<2> &indices,
                        Tensor_Options fancy_complex,
                        Weight_Options weight_function) const;

  // getter functions
  std::complex<double> get_r_p(double omega, double k) const;
  std::complex<double> get_r_s(double omega, double k) const;

  double get_za() const { return this->za; };
  double get_delta_cut() const { return this->delta_cut; };
  double get_rel_err_0() const { return this->rel_err(0); };
  double get_rel_err_1() const { return this->rel_err(1); };
  double omega_ch() const;

  // setter function
  void set_za(double za_new) { this->za = za_new; };
};
```
### Integrands
The integrands in the `GreensTensorPlate` are represented in polar coordinates such that $\mathbf{k} = (k cos(\phi), k sin(\phi))^\top$. Furthermore, the integration w.r.t $k$ was substituted with a integratio along $\kappa$. The integrand takes therefore the form
$$
\int \frac{\mathrm{d} \mathbf{k}}{2 \pi} \underline{g}(\kappa,\phi, z_a,z_a,\omega) = \frac{1}{4 \pi^2 \varepsilon_0} \int_0^\pi \mathrm{d} \phi \Big( \int_{- i \omega/c}^{-i0} + \int_0^\infty \Big) \mathrm{d} \kappa e^{-2 z_a \kappa} \Big( 1 - \frac{\cos(\phi) v \omega}{kc^2} \Big)^2 \\
 \Big\{ \kappa^2 r^p(\omega,k) \mathrm{diag} \left[ \begin{matrix} cos^2(\phi) \\ sin^2(\phi) \\ k^2/\kappa^2 \end{matrix} \right] + \frac{\omega^2}{c^2} r^s(\omega,k) \mathrm{diag} \left[ \begin{matrix} sin^2(\phi) \\ cos^2(\phi) \\ 0 \end{matrix} \right] \Big\} \;, 
\\
$$
where terms linear in $k_y$ where already neglected.
###  `# double integrand_1d_k(double phi, double omega, const vec::fixed<2> &indices,Tensor_Options fancy_complex,Weight_Options weight_function) const;`
Implements the integrand w.r.t $\phi$.

### `#double integrand_2d_k(double kappa_double, double omega, double phi,const vec::fixed<2> &indices,Tensor_Options fancy_complex,Weight_Options weight_function) const;`
Implements the integrand w.r.t $\kappa$.
                        
### `#double omega_ch();`
Returns the characteristic frequency of the macroscopic surface given by $\omega_{ch} = \frac{\delta_{cut} v}{z_a}$ with $\delta_{cut}$ being the cut of value for the $\kappa$ integration.
                        
                        

## GreensTensorPlateVacuum
Implements the sum of a [GreensTensorPlate](#GreensTensorPlate) and a [GreensTensorVacuum](#GreensTensorVacuum).
```cpp
class GreensTensorPlateVacuum : public GreensTensorPlate {
private:
  std::shared_ptr<GreensTensorVacuum> vacuum_greens_tensor;

public:
  // constructors
  GreensTensorPlateVacuum(double v, double beta, double za,
                          std::shared_ptr<ReflectionCoefficients> reflection_coefficients,
                          double delta_cut, vec::fixed<2> rel_err);
  explicit GreensTensorPlateVacuum(const std::string &input_file);

  // calculate the tensor in frequency and momentum space
  void calculate_tensor(double omega, vec::fixed<2> k,
                        cx_mat::fixed<3, 3> &GT) const override;

  // integrate over a two-dimensional k space
  void integrate_k(double omega, cx_mat::fixed<3, 3> &GT,
                   Tensor_Options fancy_complex,
                   Weight_Options weight_function) const override;

  // getters
  std::shared_ptr<GreensTensorVacuum> &get_vacuums_greens_tensor() {
    return vacuum_greens_tensor;
  };

  // setters
  void set_v(double v) override {
    this->v = v;
    this->vacuum_greens_tensor->set_v(v);
  };

};
```

For further reference on the functions see the base class [GreensTensor](#GreensTensor).

## Input file
The input file sections for the Green's tensor look like this
<!-- tabs:start -->
#### **GreensTensorVacuum**
```json
{
  "GreensTensor" : {
    "type" : "vacuum",
    "v" : ,
    "beta" : ,
    "rel_err_1" : 
  }
}
```

#### **GreensTensorPlate**
```json
{
  "GreensTensor" : {
    "type" : "plate",
    "v" : ,
    "beta" : ,
    "za" : ,
    "delta_cut" : ,
    "rel_err_0" : ,
    "rel_err_1" : 
  }
}
```
For the plate Green's tensor you also need to define the [Reflection Coefficients](api/reflection)!

#### **GreensTensorPlateVacuum**
```json
{
  "GreensTensor" : {
    "type" : "plate",
    "addvacuum" : "true",
    "v" : ,
    "beta" : ,
    "za" : ,
    "delta_cut" : ,
    "rel_err_0" : ,
    "rel_err_1" : 
  }
}
```
For the plate-and-vacuum Green's tensor you also need to define the [Reflection Coefficients](api/reflection)!
<!-- tabs:end -->

## Examples
<!-- tabs:start -->

#### ** Example : Vacuum **
We want to calculate the twofold integrated imaginary part of the vaccuum Green's tensor at frequency $\omega=3\,\mathrm{eV}$ in the static case $(v=0)$, and at $\beta=0.1\,\mathrm{eV}^{-1}$. Moreover, we need to define the accuracy of the numerical integration. Since the first integration of the vacuum Green's tensor is implemented analytically, we solely need to provide the demanded relative accuracy for the second integration $\delta_\mathrm{err}=10^{-9}$.
Since we have to define a lot of parameters, QuaCa offers a shortcut to the long task before.
We simply define all parameters in a file called `parameters.json` which looks like this
```json
{
  "GreensTensor" : {
    "type" : "vacuum",
    "v" : 0,
    "beta" : 1e-1,
    "rel_err_1" = 1e-9
  }
}
```
Now we can easily define the vacuum Green's tensor class and calculate the desired part of the integrated Green's tensor
```cpp
GreensTensorVacuum greens_tensor("parameters.json");

cx_mat::fixed<3,3> GT(fill::zeros);
greens_tensor(3.,GT, IM, UNIT);
```

#### ** Example : Plate **

We want to calculate the twofold integrated imaginary part of the scattered Green's tensor at frequency $\omega=3\,\mathrm{eV}$ in the static case $(v=0)$, at a distance of $z_a=10\,\mathrm{nm}=0.05\,\mathrm{eV}^{-1}$ (for a quick conversion between the systems of units you can use [our converter](documentation/units), and at $\beta=0.1\,\mathrm{eV}^{-1}$.

Additionally, we need to specify the geometry and material of the surface. Here, we assume a semi-infinite bulk of gold with $\omega_p=9\,\mathrm{eV}$ and $\gamma=0.1\,\mathrm{eV}$, from which we can construct the [permittivity](api/permittivity) and the [reflection coefficients](api/reflection).

 Moreover, we need to define the accuracy of the numerical integration. If we demand the relative accuracy for the second integration to be $\delta_\mathrm{rel}^1=10^{-9}$, we ought to choose a more precise relative accuracy of the integration below, as e.g. $\delta_\mathrm{rel}^0=10^{-12}$. Furthermore, in order to perform a efficient integration we provide a cut-off value for the $\kappa$ integration, since --- as shown above --- the scattered Green's tensor decays for real $\kappa$ with $\propto\exp(-2z_a\kappa)$. Since we defined the cut-off as $\kappa_\mathrm{cut}=\delta_\mathrm{cut}/(2z_a)$, a $\delta_\mathrm{cut}=30$ yields a prefactor of $\exp(-30)\approx 10^{-13}$, which should be a reasonable choice.

Since we have to define a lot of parameters, QuaCa offers a shortcut to the long task before.
We simply define all parameters in a file called `parameters.json`. Here, the permittivity and the reflection coefficients are treated as separated objects.
```json

{
    "GreensTensor": {
        "type": "plate",
        "v": 0e0,
        "za": 0.05,
        "beta": 1e-1,
        "delta_cut": 30,
        "rel_err_0": 1e-12,
        "rel_err_1": 1e-9
    },
    "ReflectionCoefficients": {
        "type": "local bulk",
    },
    "Permittivity": {
        "type": "drude",
        "gamma": 0.1,
        "omega_p": 9
    }
}
```
Now we can easily define the plate Green's tensor class and calculate the desired part of the integrated Green's tensor
```cpp
GreensTensorVacuum greens_tensor("parameters.json");

cx_mat::fixed<3,3> GT(fill::zeros);
greens_tensor.integrate_k( 3., GT, IM, UNIT);
```

#### ** Example : Plate + Vacuum **
At last we would like to create a Green's tensor with both the vacuum and the scattered contribution form the plate and again integrate the imaginary part w.r.t $\mathbf{k}$ for a frequency of $\omega = 3$. For this we simply have to add another argument to the input file we already created for the `GreensTensorPlate` which is called `addvacuum` and which we set to `true`.
```json

{
    "GreensTensor": {
        "type": "plate",
        "addvacuum": "true",
        "v": 0e0,
        "za": 0.05,
        "beta": 1e-1,
        "delta_cut": 30,
        "rel_err_0": 1e-12,
        "rel_err_1": 1e-9
    },
    "ReflectionCoefficients": {
        "type": "local bulk",
    },
    "Permittivity": {
        "type": "drude",
        "gamma": 0.1,
        "omega_p": 9
    }
}
```
Storing the above arguments in a file called `parameter.json` we can now easily compute the twofold integration w.r.t. $\mathbf{k}$ for the imaginary part of the full Green's tensor. We just have to create an instance of `GreensTensorPlateVacuum` instead of `GreensTensorPlate`.
```cpp
  GreensTensorPlateVacuum greens_tensor("parameter.json");
  
  greens_tensor.integrate_k(3., GT, IM, UNIT);
```


<!-- tabs:end -->
