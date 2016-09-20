
// This is an example in README.md.

#include <iostream>
#include "units.hpp"
using namespace num;
using namespace num::u;
int main() {
   using time = num::time;  // in case of conflict with C library's time()
   length x0 = 2 * m;       // meters
   speed  v  = 5 * m / s;   // meters per second
   time   t  = 3 * min;     // minutes
   std::cout << x0 + v * t << std::endl;
   std::cout << (x0 + v * t) / yd << " yards" << std::endl;
   return 0;
}

