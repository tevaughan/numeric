
/// \file   interp_1.cpp
/// \brief  Example of how to use interpolant.

#include <cmath>
#include <iostream>

#include "interpolant.hpp"
#include "integral.hpp"

using namespace num;
using namespace std;

/// Gaussian with zero mean and unitary standard deviation.
double g(double x)
{
   return exp(-0.5 * x * x) / sqrt(2.0 * M_PI);
}

/// Main function for interp_1.cpp.
int main()
{
   interpolantd i;
   integral_rk(g, -5.0, 5.0, 1.0E-04, 2, &i);
   cout << i.points();
   return 0;
}

