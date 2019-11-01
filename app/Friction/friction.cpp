#include "Quaca.h"
#include <iostream>

int main(int argc, char *argv[])
{
  Polarizability *pol = new Polarizability();
  double test = pol->muR(5.0);

  std::cout << "Hello adorable mammal!" << std::endl;
  std::cout << test << std::endl;
  return 0;
};
