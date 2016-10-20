
/// \file   interp_1.cpp
/// \brief  Example of how to use interpolant.

#include <cmath>
#include <fstream>
#include <iostream>

#include "rk.hpp"

using namespace num;
using namespace std;

/// Gaussian with zero mean and unitary standard deviation.
double g(double x) { return exp(-0.5 * x * x) / sqrt(2.0 * M_PI); }

/// Main function for interp_1.cpp.
int main()
{
   double const x1 = -5.0;
   double const x2 = +5.0;
   rk_quadd const i(g, x1, x2, 1.0E-04, 2, true);
   cout << "fnc_points:\n" << i.intermed_fnc() << endl;
   cout << "int_points:\n" << i.intermed_int() << endl;
   ofstream os("interp_2.dat");
   ofstream osa("interp_2a.dat");
   ofstream osb("interp_2b.dat");
   os << i.intermed_fnc() << endl;
   unsigned const n = 1000;
   double const dx = (x2 - x1)/(n - 1);
   auto const interp = i.make_fnc_interp();
   for (unsigned j = 0; j < n; ++j) {
      double const x = x1 + dx * j;
      osa << x << " " << interp(x) << endl;
      osb << x << " " << g(x) << endl;
   }
   return 0;
}

