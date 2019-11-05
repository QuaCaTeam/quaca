#include "catch.hpp"
#include "Quaca.h"
#include <iostream>

double f( double x, void *p)
{
  return x*x;
};

  TEST_CASE("Integration test")
{
double test=cquad(f,0.0,1.0,1E-5,0);
std::cout << test << std::endl;
};
