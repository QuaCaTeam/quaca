# Implement a new feature
In this tutorial, we will add a new feature to QuaCa.
More specifically, we define a new permittivity function, i.e. the Lorentz model of the form
$$
\varepsilon(\omega) = \varepsilon_{\infty} - \frac{\alpha_0 \omega_0^2}{\omega_0^2 - \omega^2 - i \omega \gamma}
$$,
where $\varepsilon_{\infty}$, $\alpha_0$, $\omega_0^2$, and $\gamma$ are constants.

To this end, we need to:

1. Write tests for the permittivity
2. Implement the permittivity
3. Check if the tests succeed
4. (Document the code)

Let us first look at how a generic permittivity looks like in QuaCa.
Move to the directory `src/Permittivity/` and open the file `Permittivity.h`.
You will see the following code
```cpp
#ifndef PERMITTIVITY_H
#define PERMITTIVITY_H

#include <complex>

//! An abstract permittivity class
class Permittivity {
public:
  // calculate the permittivity
  virtual std::complex<double> calculate(double omega) = 0;

  // calculate the permittivity times omega
  virtual std::complex<double> calculate_times_omega(double omega) = 0;

  // print info
  virtual void print_info(std::ostream &stream) const =0;
};

#endif // PERMITTIVITY_H
```
This is a so-called abstract class.
It is abstract in the sense that all functions defined within the class are preceded by the prefix `virtual` and appended by the expression `=0`.
No specific model is needed at the moment as the syntax allows for a certain level of generality in defining the functions.
All permittivities are children of this class, so they inherit from it and they have to concretely implement the three functions `calculate (double omega)`, `calculate_times_omega (double omega)` and `print_info(std:ostream &stream)`.

Every class also needs a so-called constructor that defines the parameters of the class (in our case, e.g., $\varepsilon_{\infty}$).

Let us now start by writing tests for our new permittivity.

## 1. Write tests for the permittivity
In the QuaCa project, we are practicing so-called [test driven development](dev/testing).
Before we implement a new feature, we will write (physically motivated) tests that this feature has to pass.
Given a permittivity, for instance, we might want to check whether it obeys the crossing relation
$$
\varepsilon(\omega) = \varepsilon^{*}(-\omega)
$$.

We now describe the general procedure of implementing tests step by step using the constructors of the class.
We will turn to the concrete example of testing the crossing relation later.
Change to the directory `test/UnitTests/Permittivity` and create a new file called `test_PermittivityLorentzOhmic_unit.cpp` (PermittivityLorentzOhmic will be the name of the class we create later).
Fill this file with the following lines
```cpp
#include "Quaca.h"
#include "catch.hpp"

TEST_CASE("LorentzOhmic permittivity constructors work as expected",
          "[PermittivityLorentzOhmic]") {
};
```
We simply test the constructors of our class.
In the last two tutorials, you have seen that QuaCa can process parameter settings in two ways: Either you give it an input file or you define your own main function.
Therefore, we want to implement two constructors to check both of these features.
One takes all the arguments directly
```cpp
PermittivityLorentzOhmic(double eps_inf, double alpha_0, double omega_0, double gamma);
```
And the other one takes a path to an input file and reads the parameters from this file
```cpp
PermittivityLorentzOhmic(std::string input_file);
```
We now test the functionality of these two constructors.
Create a section for the direct constructor and define our permittivity
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
Now check whether the given parameters were processed correctly
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
    REQUIRE(perm.get_alpha_0() == alpha_0);
    REQUIRE(perm.get_omega_0() == omega_0);
    REQUIRE(perm.get_gamma() == gamma);
  };
};
```
When we have implemented the constructors and the getter functions correctly, this test will succeed and we can be sure that our new feature is constructed properly.

Let us now advance to a physically relevant example introduced above and check whether our permittivity obeys the crossing relation.
We add a new test case below the constructor test:
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
    REQUIRE(perm.get_alpha_0() == alpha_0);
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

  auto omega = GENERATE(-1.2,0.1);
  REQUIRE(perm.calculate(omega) == std::conj(perm.calculate(-omega)));
};
```
This will test the crossing relation for the set of frequencies defined between the brackets after the `GENERATE` command.
In practice, be sure to include tests in all the relevant regimes of the function.
In our case, e.g., we want to ensure that the crossing relation is fulfilled for positive and negative frequencies.
As a last step, we have to add this test to the QuaCa functionalities.
Change to the `test/UnitTests` directory and open `CMakeLists.txt`.
You will see a number of files that are added to the target `test_sources`.
Add our test to the end of this list
```cmake
# add sources to test
set(test_sources
        ...
        Permittivity/test_PermittivityDrude_unit.cpp
        Permittivity/test_PermittivityLorentz_unit.cpp
        Permittivity/test_PermittivityLorentzOhmic_unit.cpp # This is the line you have to add
        ...
        )
```
The next time we compile, CMake will know about this test and weave it into our test executable.
With the tests ready, let us now actually implement the new feature.

## 2. Implement the permittivity
Change to the `src/Permittivity/` directory and create a new file called `PermittivityLorentzOhmic.h`.
Add the following code
```cpp
#ifndef PERMITTIVITYLORENTZOHMIC_H
#define PERMITTIVITYLORENTZOHMIC_H

#include "Permittivity.h"
#include <complex>

//! A LorentzOhmic model permittivity
class PermittivityLorentzOhmic : public Permittivity {
private:
  double eps_inf;
  double alpha_0;
  double omega_0;
  double gamma;

public:
  // constructors
  PermittivityLorentzOhmic(double eps_inf, double alpha_zero, double omega_0,
                           double gamma);
  PermittivityLorentzOhmic(const std::string &input_file);

  // calculate the permittivity
  std::complex<double> calculate(double omega) const override;
  // calculate the permittivy times omega
  std::complex<double> calculate_times_omega(double omega) const override;

  // getter methods
  double get_eps_inf() const { return this->eps_inf; };
  double get_alpha_0() const { return this->alpha_0; };
  double get_omega_0() const { return this->omega_0; };
  double get_gamma() const { return this->gamma; };

  // print all physical quantities
  void print_info(std::ostream &stream) const override;
};
#endif
```
The header file defines all functions that a class contains. We have to define the actual implementation of the physical model in the additional file `src/Permittivity/PermittivityLorentzOhmic.cpp`:

```cpp
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "PermittivityLorentzOhmic.h"

PermittivityLorentzOhmic::PermittivityLorentzOhmic(double eps_inf,
                                                   double alpha_0,
                                                   double omega_0, double gamma)
    : eps_inf(eps_inf), alpha_0(alpha_0), omega_0(omega_0), gamma(gamma) {
}

PermittivityLorentzOhmic::PermittivityLorentzOhmic(
    const std::string &input_file) {
  // Load the json file in ptree
  pt::ptree root;
  pt::read_json(input_file, root);

  // check if type is right
  std::string type = root.get<std::string>("Permittivity.type");
  assert(type == "lorentz ohmic");

  // read parameters
  this->eps_inf = root.get<double>("Permittivity.eps_inf");
  this->alpha_0 = root.get<double>("Permittivity.alpha_0");
  this->omega_0 = root.get<double>("Permittivity.omega_0");
  this->gamma = root.get<double>("Permittivity.gamma");
}

std::complex<double> PermittivityLorentzOhmic::calculate(double omega) const {
  // dummies for result and complex unit
  std::complex<double> result;
  std::complex<double> I(0.0, 1.0);

  // calculate the result
  result =
      eps_inf - alpha_0 * omega_0 * omega_0 /
                    (omega_0 * omega_0 - omega * omega - I * gamma * omega);

  return result;
}

std::complex<double> PermittivityLorentzOhmic::calculate_times_omega(double omega) const {
  // dummies for result and complex unit
  std::complex<double> result;
  std::complex<double> I(0.0, 1.0);

  // calculate the result
  result =
      eps_inf - alpha_0 * omega_0 * omega_0 * omega /
                    (omega_0 * omega_0 - omega * omega - I * gamma * omega);

  return result;
}

void PermittivityLorentzOhmic::print_info(std::ostream &stream) const {
  stream << "# PermittivityLorentzOhmic\n#\n"
         << "# eps_inf = " << eps_inf << "\n"
         << "# alpha_0 = " << alpha_0 << "\n"
         << "# omega_0 = " << omega_0 << "\n"
         << "# gamma = " << gamma << "\n";
}
```

Again, we have to tell QuaCa to compile this target.
Change to the `src/` directory and open the `CMakeLists.txt`.
Add the source file `PermittivityLorentzOhmic.cpp` to the other source files
```cmake
# add QuaCa library
set(quaca_sources
        ...
        Permittivity/PermittivityDrude.cpp
        Permittivity/PermittivityFactory.cpp
        Permittivity/PermittivityLorentz.cpp
        Permittivity/PermittivityLorentzOhmic.cpp # This is the line you have to add
        ...
        )
```
You also need to include the header file `PermittivityLorentzOhmic.h` in the public header file `QuaCa.h`, located in the `include/` directory
```cpp
#ifndef QUACA_H
#define QUACA_H

...

#include "../src/Permittivity/Permittivity.h"
#include "../src/Permittivity/PermittivityDrude.h"
#include "../src/Permittivity/PermittivityFactory.h"
#include "../src/Permittivity/PermittivityLorentz.h"
#include "../src/Permittivity/PermittivityLorentzOhmic.h" // This is the line you have to add

...

#endif //QUACA_H
```

# 3. Check if the tests succeed
Build the project by typing in the `build/` directory
```bash
build/ $ cmake ..
build/ $ make
```
and run the unit tests by typing
```bash
build/ $ ./../bin/test_quaca_unit
```

If you followed the above steps, there should be no error, showing us that our tests indeed succeeded.

# (4. Document the code)
For resons of brevity, we did not include comments in the above code examples.
You should, however, include documentation on what you have coded in the form of comments inside the header files
as well as in the source files.
Ideally, you should also write a corresponding `.md` file in the documentation folder `docs/`.
