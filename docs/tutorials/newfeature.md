!> TODO: finish writing and test this tutorial

# Implement a new feature
For now we have played around with the main functions and executables.
In this tutorial we will add a feature to QuaCa, calculate the quantum friction with the new feature and see what changes.

The new feature we are going to implement is a new permittivity.
More specifically we want to implement a simple Lorentz model of the form
$$
\varepsilon(\omega) = \varepsilon_{\infty} - \frac{\alpha_0 \omega_0^2}{\omega_0^2 - \omega^2 - i \omega \gamma}.
$$

For this we need to:

1. Write tests for the permittivity
2. Implement the permittivity
3. Check if the tests succeed
4. Calculate the friction with the new permittivity


Let us first look at how a generic permittivity looks like in QuaCa.
For this head over to the directory `src/Permittivity/` and open the file `Permittivity.h`.
You will see the following code
```cpp
#ifndef PERMITTIVITY_H
#define PERMITTIVITY_H

#include <complex>

//! An abstract permittivity class
class Permittivity {
public:
  // calculate the permittivity
  virtual std::complex<double> epsilon(double omega) = 0;

  // calculate the permittivity times omega
  virtual std::complex<double> epsilon_omega(double omega) = 0;
};

#endif // PERMITTIVITY_H
```
This is a so called abstract class, it is abstract because in front of all functions that are defined within the class you find the word `virtual` and at the end you find `=0`.
This means that we do not know that the function in general actually does, which makes sense because we simply do not know what the specific permittivity is for now (i.e. how we model the permittivity).
All permittivities are children of this class, so they inherit from it and they have to implement the two functions `epsilon (double omega)` and `epsilon_omega (double omega)`.

Every class also needs a so called constructors that basically defines the parameters of the class (in our case e.g. the blabla $\varepsilon_{\infty}$).

Let us now start by writing tests for our new permittivity.

## 1. Write tests for the permittivity
In the QuaCa project we are practicing so called [test driven development](testing).
This means that before we implement a new feature we will write tests that this feature has to fulfill.
Seeing that we want to implement a permittivity, we might want to check whether we construct it correctly and whether it really obeys the crossing relation
$$
\varepsilon(\omega) = \varepsilon^{*}(-\omega).
$$

Let us start by heading over to the directory `test/UnitTests` and creating a new file called `test_PermittivityLorentzOhmic_unit.cpp` (PermittivityLorentzOhmic will be the name of the class we create later).
Fill this file with the following code
```cpp
#include "Quaca.h"
#include "catch.hpp"

TEST_CASE("LorentzOhmic permittivity constructors work as expected",
          "[PermittivityLorentzOhmic]") {
};
```
In this test case we want to test the constructors of our class.
In the last two tutorials you have seen that QuaCa works in two ways: either you give it an input file or you define your own main function.
Therefore we want to implement two constructors that look like this: one takes all the arguments directly
```cpp
PermittivityLorentzOhmic(double eps_inf, double alpha_0, double omega_0, double gamma);
```
And the other one takes a path to an input file and reads the parameters from this file
```cpp
PermittivityLorentzOhmic(std::string input_file);
```
How would we tests if this actually works?
Well, we simply define the quantities and check if they are the want we defined.
Let us create a section for the direct constructor and define our permittivity
```cpp
#include "Quaca.h"
#include "catch.hpp"

TEST_CASE("LorentzOhmic permittivity constructors work as expected",
          "[PermittivityLorentzOhmic]") {
  SECTION("Direct constructor") {
    double eps_inf = 1.4;
    double alpha_0 = 6e-9;
    double omega_0 = 3.4;
    double gamma = 2.8;

    PermittivityLorentzOhmic perm(eps_inf, alpha_0, omega_0, gamma);
  };
};
```
Now check whether the given parameters were given correctly
```cpp
#include "Quaca.h"
#include "catch.hpp"

TEST_CASE("LorentzOhmic permittivity constructors work as expected",
          "[PermittivityLorentzOhmic]") {
  SECTION("Direct constructor") {
    double eps_inf = 1.4;
    double alpha_0 = 6e-9;
    double omega_0 = 3.4;
    double gamma = 2.8;

    PermittivityLorentzOhmic perm(eps_inf, alpha_0, omega_0, gamma);

    REQUIRE(perm.get_eps_inf() == eps_inf);
    REQUIRE(perm.get_alpha_zero() == alpha_zero);
    REQUIRE(perm.get_omega_0() == omega_0);
    REQUIRE(perm.get_gamma() == gamma);
  };
};
```
And you are done with this test!
When we have implemented our permittivity correctly, i.e. we have implemented the constructors and the getter functions correctly this test will succeed and we can be sure that our new feature is constructed correctly.

Let us now check whether our permittivity obeys the crossing relation in a new test case below the constructor test
```cpp
#include "Quaca.h"
#include "catch.hpp"

TEST_CASE("LorentzOhmic permittivity constructors work as expected",
          "[PermittivityLorentzOhmic]") {
  SECTION("Direct constructor") {
    double eps_inf = 1.4;
    double alpha_0 = 6e-9;
    double omega_0 = 3.4;
    double gamma = 2.8;

    PermittivityLorentzOhmic perm(eps_inf, alpha_0, omega_0, gamma);

    REQUIRE(perm.get_eps_inf() == eps_inf);
    REQUIRE(perm.get_alpha_zero() == alpha_zero);
    REQUIRE(perm.get_omega_0() == omega_0);
    REQUIRE(perm.get_gamma() == gamma);
  };
};

TEST_CASE("LorentzOhmic permittivity obeys crossing relation",
          "[PermittivityLorentzOhmic]") {
  double eps_inf = 1.4;
  double alpha_0 = 6e-9;
  double omega_0 = 3.4;
  double gamma = 2.8;

  PermittivityLorentzOhmic perm(eps_inf, alpha_0, omega_0, gamma);

  auto omega = GENERATE(take(10, random(-150.4, 150.4)));
  REQUIRE(perm.epsilon(omega) == std::conj(perm.epsilon(-omega)));
};
```
And done!
This test will generate ten randomly drawn frequencies and check whether the crossing relation is fulfilled for them.
With the tests ready to run let us now actually implement the new feature.

## 2. Implement the permittivity
Head over to the `src/Permittivity` directory and create a new file called `PermittivityLorentzOhmic.h`.
Fill it with the following code
```cpp
#ifndef PERMITTIVITYLORENTZOHMIC_H
#define PERMITTIVITYLORENTZOHMIC_H

#include "Permittivity.h"
#include <complex>

//! A LorentzOhmic model permittivity
class PermittivityLorentzOhmic : public Permittivity {
private:
  double eps_inf;
  double alpha_zero;
  double omega_0;
  double gamma;

public:
  // constructors
  PermittivityLorentzOhmic(double eps_inf, double alpha_zero, double omega_0, double gamma);
  PermittivityLorentzOhmic(std::string input_file);

  // calculate the permittivity
  std::complex<double> epsilon(double omega);

  // Returns the numerical value of the permittivity scaled by omega.
  std::complex<double> epsilon_omega(double omega);

  // getter methods
  double get_eps_inf() { return this->eps_inf; };
  double get_alpha_0() { return this->alpha_0; };
  double get_omega_0() { return this->omega_0; };
  double get_gamma() { return this->gamma; };
};

#endif // PERMITTIVITYLORENTZOHMIC_H
```
The header file defines all functions that a class contains.
