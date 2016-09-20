
/// \file   units_1.cpp
/// \brief  Example of how to use units.

#include <iostream>
#include "units.hpp"

using namespace num;
using namespace num::u;
using namespace std;

/// Main function for units_1.cpp.
int main()
{
   // 'ft' stands for one foot.
   length const ell = 2.1 * ft;        // length of pendulum
   // 'gn' stands for one gee.
   frequency const f = sqrt(gn / ell); // frequency of small oscillations
   cout << f << endl;
   return 0;
}

