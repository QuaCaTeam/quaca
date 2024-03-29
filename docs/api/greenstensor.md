# GreensTensor {docsify-ignore-all}
This is an abstract class that defines a Green's tensor.
A specific kind of Green's tensor will be a child of this class.

## Member functions

### `GreensTensor(std::string input_file)`
Description: Input file constructor of the class.
* Input parameters: 
    - `std::string input_file`: Name of the input file
* Return value: `void`

### `GreensTensor(double v, double beta)`
Direct constructor of the class.
* Input parameters:
    - `double v`: velocity of the atom
    - `double beta`: inverse temperature of the system
* Return value: `void`:

### `virtual void calculate_tensor(double omega, vec::fixed<2> k,cx_mat::fixed<3, 3> &GT) const = 0;`
Calculates the Green's tensor, i.e. $\underline{G}(\mathbf{k},z,z', \omega + k_x v)$, where $x$ is the axis along the motion with velocity $v$, and puts the result into the matrix `GT`. Here, the $x$ and $y$ coordinate are calculated in reciprocal space $(x-x',y-y') \to (k_x,k_y)$, while the $z$ or $z'$ component remain in the space domain. The $z$ components are defined in the respective child class. For more physical details, see [this publication](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.94.042114). 
* Input parameters:
    - `double omega`: Frequency, at which the Green's tensor is evaluated.
    - `vec:fixed<2> k`: momentum vector $\mathbf{k} = (k_x, k_y)$.
    - `cx_mat:fixed<3, 3> &GT`: 3x3 matrix, where the resulting Green's tensor is stored.
* Return value: `void`

###  `virtual void integrate_k(double omega, cx_mat::fixed<3, 3> &GT, Tensor_Options fancy_complex, Weight_Options weight_function) const = 0;`
Compute the integral over the momentum space of the Green's tensor i.e. $\int_{-\infty}^{\infty} \frac{\mathrm{d}^2 \mathbf{k}}{(2 \pi)^2} \; f(\mathbf{k},\omega) \underline{G}(\mathbf{k},z,z', \omega + k_x v)$, where $f(\mathbf{k},\omega)$ is some weight function.  The result is stored in `GT`.

* Input parameters:
    - `double omega`: Frequency, at which the Green's tensor is evaluated.
    - `cx_mat:fixed<3, 3> &GT`: 3x3 matrix, where the resulting Green's tensor is stored.
    - `Tensor_Options fancy_complex`: options arel isted below
    - `Weight_Options weight_function`: options are listed below
* Return value: `void`

The different options for the two variables `fancy_complex` and `weight_function` are

| fancy_complex                | Math       |
|---------------------|------------|
| `RE`           | $\underline{G}_\Re = (\underline{G}+\underline{G}^\dagger)/2$ |
| `IM`           | $\underline{G}_\Im = (\underline{G}-\underline{G}^\dagger)/(2\mathrm{i})$ |

and

| weight_function                | Math       |
|---------------------|------------|
| `UNIT` | $\times \mathbb{1}$ | 
| `KV`           | $\times \, k_x $ |
| `TEMP`           | $\times\left[1-\exp(-\hbar\beta(\omega+k_xv))\right]^{-1} $ |
| `KV_TEMP`           | $\times k_x\left[1-\exp(-\hbar\beta(\omega+k_xv))\right]^{-1} $ |
| `NON_LTE`           | $\times \left( \left[1-\exp(-\hbar\beta(\omega+k_xv))\right] - \left[1-\exp(-\hbar\beta\omega) \right]  \right) $ |
| `NON_LTE`           | $\times k_x \left( \left[1-\exp(-\hbar\beta(\omega+k_xv))\right] - \left[1-\exp(-\hbar\beta\omega) \right]  \right) $ |


### `# virtual double omega_ch() const = 0;`
Calculates and returns a characteristic frequency $\omega_\mathrm{ch}$ of the respective Green's tensor
* Input parameters: `void`
* Return value:
    - `double`: value of the characteristic frequency

### `# double get_...`
Are getter functions that return the respective quantity (here, for example, `v` or `beta`)
* Input parameters: `void`
* Return value: `mixed`

### `# virtual void set_...`
Are setter functions which enable to set the respective attribute of the class (here `v`) to a desired value.
* Input parameters: `mixed`
* Return value: `void`

# GreensTensorVacuum
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
where $\omega$ is the frequency, $\mathbf{k}=(k_x,\,k_y)^\intercal$ (with $|\mathbf{k}|=k$) is the corresponding wavevector to the two-dimensional $(x-x')\wedge (y-y')$ plane, and $z-z'$ is the remaining (not Fourier-transformed) spatial coordinate. The $z- z'\to 0$ indicates, that the Green's tensor is evaluated twice at the same point with respect to the $z$ coordinate. The introduced Heaviside function $\theta(\omega^2/c^2-k^2)$ ensures, that solely imaginary values are considered. The real part of the vacuum Green's tensor is implicitly contained in the $\omega_a$. A clear derivation of the vacuum Green's tensor can be found for instance in [this paper](https://www.mdpi.com/2076-3417/7/11/1158).

## Member functions

### `double integrand_k(double kv, double omega, const uvec::fixed<2> &indices,Tensor_Options fancy_complex,Weight_Options weight_function) const;`
Implements the integrand in momentum space along the direction of motion of the microscopic particle.
* Input parameters:
    - `double kv`: momentum in the direction of motion
    - `double omega`: frequency
    - `const uvec::fixed<2> &indices`: Indices for entries of the Green's tensor, that should be computed.
    - `Tensor_Options fancy_complex`: See [GreensTensor](#GreensTensor) for options.
    - `Weight_Options weight_function`: See [GreensTensor](#GreensTensor) for options.
* Return value:
    - `double`: value of the integrand

# GreensTensorPlate
Implements the scattered part of the Green's tensor of a flat surface beneath vacuum and is given by
$$
  \underline{g}(\mathbf{k},z_a, z_a,\omega) =
  \frac{\omega^2}{c^2}
  \big(
  r^s(\mathbf{k},\omega) \, \mathbf{ e}_s \otimes {\mathbf{ e}_s}
  +
  r^p(\mathbf{k},\omega) \,\cdot \mathbf{ e}_{+p} \otimes {\mathbf{ e}_{-p}}
\big)\frac{e^{-\kappa (z_a+z_a)}}{2\epsilon_0\kappa}
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

## Integrands
The integrands in the `GreensTensorPlate` are represented in polar coordinates such that $\mathbf{k} = (k \cos\phi, k \sin\phi)^\top$. Furthermore, the integration with respect to $k$ was substituted with a integration along $\kappa$. The integrand takes therefore the form
$$
\int \frac{\mathrm{d}^2 \mathbf{k}}{(2 \pi)^2} \underline{g}(\kappa,\phi, z_a,z_a,\omega) = \frac{1}{4 \pi^2 \varepsilon_0} \int_0^\pi \mathrm{d} \phi \Big( \int_{- \mathrm{i} \omega/c}^{-\mathrm{i}0} + \int_0^\infty \Big) \mathrm{d} \kappa e^{-2 z_a \kappa} \Big( 1 - \frac{v \omega\cos\phi }{kc^2} \Big)^2 \\
 \Big\{ \kappa^2 r^p(\omega,k) \mathrm{diag} \left[ \begin{matrix} \cos^2\phi \\ \sin^2\phi \\ k^2/\kappa^2 \end{matrix} \right] + \frac{\omega^2}{c^2} r^s(\omega,k) \mathrm{diag} \left[ \begin{matrix} \sin^2\phi \\ \cos^2\phi \\ 0 \end{matrix} \right] \Big\} \;, 
\\
$$
where terms linear in $k_y$ where already neglected.

## Member functions
###  `# double integrand_1d_k(double phi, double omega, const uvec::fixed<2> &indices,Tensor_Options fancy_complex,Weight_Options weight_function) const;`
Implements the integrand with respect to $\phi$.
* Input parameters:
    - `double phi`: angle of the polar coordinates
    - `double omega`: frequency of the integrand
    - `const uvec::fixed<2> &indices`: indices of the component of the Green's tensor, which should be computed
    - `Tensor_Options fancy_complex`: see [GreensTensor](#GreensTensor) for options
    - `Weight_Options weight_function`: see [GreensTensor](#GreensTensor) for options
* Return value:
    - `double`: value of the integrand


### `#double integrand_2d_k(double kappa_double, double omega, double phi,const uvec::fixed<2> &indices,Tensor_Options fancy_complex,Weight_Options weight_function) const;`
Implements the integrand with respect to $\kappa$.
* Input parameters:
    - `double kappa_double`: value of $\kappa$
    - `double omega`: frequency of the integrand
    - `double phi`: angle of the polar coordinates
    - `const uvec::fixed<2> &indices`: indices of the component of the Green's tensor, which should be computed
    - `Tensor_Options fancy_complex`: see [GreensTensor](#GreensTensor) for options
    - `Weight_Options weight_function`: see [GreensTensor](#GreensTensor) for options
* Return value:
    - `double`: value of the integrand


### `#double omega_ch();`
Returns the characteristic frequency of the macroscopic surface given by $\omega_\mathrm{ch} = \frac{\delta_\mathrm{cut} v}{z_a}$ with $\delta_\mathrm{cut}$ being the cut of value for the $\kappa$ integration.
* Input parameters: `void`
* Return value: 
    - `double`: value of the characteristic freqency $\omega_\mathrm{ch}$

# GreensTensorPlateVacuum
Implements the sum of a [GreensTensorPlate](#GreensTensorPlate) and a [GreensTensorVacuum](#GreensTensorVacuum).
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
We want to calculate the twofold integrated imaginary part of the vacuum Green's tensor at frequency $\omega=3\,\mathrm{eV}$ in the static case $(v=0)$, and at $\beta=0.1\,\mathrm{eV}^{-1}$. Moreover, we need to define the accuracy of the numerical integration. Since the first integration of the vacuum Green's tensor is implemented analytically, we solely need to provide the demanded relative accuracy for the second integration $\delta_\mathrm{err}=10^{-9}$.
We define all needed parameters, define the vacuum Green's tensor class and calculate the desired part of the integrated Green's tensor
```cpp
// Set parameters
double v = 0;
double beta = 1e-1;
double relerr = 1e-9;
double omega = 3e0;

// Create Green's tensor class
GreensTensorVacuum greens_tensor(v, beta, relerr);

// Calculate desired Green's tensor
cx_mat::fixed<3,3> GT(fill::zeros);
greens_tensor.integrate_k(omega, GT, IM, UNIT);
```

#### ** Example (.json) : Vacuum **
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
// Create Green's tensor class
GreensTensorVacuum greens_tensor("parameters.json");

// Calculate desired Green's tensor
cx_mat::fixed<3,3> GT(fill::zeros);
greens_tensor.integrate_k(omega, GT, IM, UNIT);
```

#### ** Example : Plate **

We want to calculate the twofold integrated imaginary part of the scattered Green's tensor at frequency $\omega=3\,\mathrm{eV}$ in the static case $(v=0)$, at a distance of $z_a=10\,\mathrm{nm}=0.05\,\mathrm{eV}^{-1}$ (for a quick conversion between the systems of units you can use [our converter](documentation/units), and at $\beta=0.1\,\mathrm{eV}^{-1}$.

Additionally, we need to specify the geometry and material of the surface. Here, we assume a semi-infinite bulk of gold with $\omega_p=9\,\mathrm{eV}$ and $\gamma=0.1\,\mathrm{eV}$, from which we can construct the [permittivity](api/permittivity) and the [reflection coefficients](api/reflection).

 Moreover, we need to define the accuracy of the numerical integration. If we demand the relative accuracy for the second integration to be $\delta_\mathrm{rel}^1=10^{-9}$, we ought to choose a more precise relative accuracy of the integration below, as e.g. $\delta_\mathrm{rel}^0=10^{-12}$. Furthermore, in order to perform a efficient integration we provide a cut-off value for the $\kappa$ integration, since --- as shown above --- the scattered Green's tensor decays for real $\kappa$ with $\propto\exp(-2z_a\kappa)$. Since we defined the cut-off as $\kappa_\mathrm{cut}=\delta_\mathrm{cut}/(2z_a)$, a $\delta_\mathrm{cut}=30$ yields a prefactor of $\exp(-30)\approx 10^{-13}$, which should be a reasonable choice.

First, we define all parameters and construct the permittivity and the reflection coefficients.
```cpp
// Define parameters and generate permittivity
// and reflection coefficient class
double gamma = 0.1;
double omega_p = 9;
auto perm = std::make_shared<PermittivityDrude>(omega_p, gamma);
auto refl = std::make_shared<ReflectionCoefficientsLocBulk>(perm);       
``` 
Now we can use the remaining parameters and the generated classes to create a Green's tensor class of the plate configuration and calculate the desired part
```cpp
// Set parameters
double v = 0e0;
double za = 0.05;
double omega = 3e0;
double beta= 1e-1;
double delta_cut 30;
vec::fixed <2> rel_err = {1e-12, 1e-9};

GreensTensorPlate greens_tensor(v, beta, za, refl, delta_cut, rel_err);

cx_mat::fixed<3,3> GT(fill::zeros);
greens_tensor.integrate_k( omega, GT, IM, UNIT);
```

#### ** Example (.json): Plate **

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
At last we would like to create a Green's tensor with both the vacuum and the scattered contribution form the plate and again integrate the imaginary part with respect to $\mathbf{k}$ for a frequency of $\omega = 3\,\mathrm{eV}$. For this we simply repeat the same initialization as for the `GreensTensorPlate` and solely alter the child class. 
First, we define all parameters and construct the permittivity and the reflection coefficients.
```cpp
// Define parameters and generate permittivity
// and reflection coefficient class
double gamma = 0.1;
double omega_p = 9;
auto perm = std::make_shared<PermittivityDrude>(omega_p, gamma);
auto refl = std::make_shared<ReflectionCoefficientsLocBulk>(perm);       
``` 
Now we can use the remaining parameters and the generated classes to create a Green's tensor class of the plate configuration and calculate the desired part
```cpp
// Set parameters
double v = 0e0;
double za = 0.05;
double omega = 3e0;
double beta= 1e-1;
double delta_cut 30;
vec::fixed <2> rel_err = {1e-12, 1e-9};

GreensTensorPlateVacuum greens_tensor(v, beta, za, refl, delta_cut, rel_err);

cx_mat::fixed<3,3> GT(fill::zeros);
greens_tensor.integrate_k( omega, GT, IM, UNIT);
```


#### ** Example (.json) : Plate + Vacuum **
At last we would like to create a Green's tensor with both the vacuum and the scattered contribution form the plate and again integrate the imaginary part with respect to$\mathbf{k}$ for a frequency of $\omega = 3\,\mathrm{eV}$. For this we simply have to add another argument to the input file we already created for the `GreensTensorPlate` which is called `addvacuum` and which we set to `true`.
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
Storing the above arguments in a file called `parameter.json` we can now easily compute the twofold integration with respect to $\mathbf{k}$ for the imaginary part of the full Green's tensor. We just have to create an instance of `GreensTensorPlateVacuum` instead of `GreensTensorPlate`.
```cpp
  GreensTensorPlateVacuum greens_tensor("parameter.json");
  
  greens_tensor.integrate_k(3., GT, IM, UNIT);
```


<!-- tabs:end -->
