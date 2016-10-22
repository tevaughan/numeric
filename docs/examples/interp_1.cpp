
/// \file   interp_1.cpp
/// \brief  Example of how to use interpolant.

#include <cmath>
#include <iostream>

#include "interpolant.hpp"

using namespace num;
using namespace std;

/// Gaussian with zero mean and unitary standard deviation.
double g(double x) { return exp(-0.5 * x * x) / sqrt(2.0 * M_PI); }

/// Main function for interp_1.cpp.
int main()
{
   auto const i = make_linear_interp(g, -5.0, 5.0, 1.0E-04);
   for (auto p = i.dat().begin(); p != i.dat().end(); ++p) {
      auto const x = p->a - 0.5 * p->da;
      cout << x << " " << i(x) << endl;
   }
   return 0;
}

