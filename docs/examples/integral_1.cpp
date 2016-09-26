
// This example is in README.md.

#include <cmath>    // for exp()
#include <iostream> // for cout, endl

#include "integral.hpp"
#include "units.hpp"

using namespace num;
using namespace num::u;

using prob_density = decltype(1.0 / length());

length sigma = 3 * cm;

prob_density gaussian(length x)
{
   double const z = x / sigma;
   return exp(-0.5 * z * z) / sigma / sqrt(2.0 * M_PI);
}

int main()
{
   double const prob = integral(gaussian, -3 * cm, 3 * cm);
   std::cout << prob << std::endl; // 0.682689, about 68% for one sigma
   return 0;
}

