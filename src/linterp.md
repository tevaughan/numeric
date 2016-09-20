
# Linear Interpolation {#linterp}

A linear interpolant can be initialized directly from a file with two
space-delimited columns.  Then it can be scaled, multiplied by another
interpolant with no loss of information, integrated, etc.

```cpp
#include <iostream>
#include <interpolant.hpp>
#include <units.hpp>

using namespace num;
using namespace std;

// Read in first two columns of space-delimited file.
// First column is interpreted as a double-precision length in nanometers.
// Second column is interpreted as a double-precision time in seconds.
// Explicit namespace scoping for num::time might be required to avoid
// conflict with C library's time() function.
interpolant<length, num::time> i("file.txt", u::nm, u::s);

// Read in first two columns of space-delimited file.
// Each column is interpreted as a double.
interpolantd j("file.txt");

int main()
{
   cout << i(500 * u::nm) << endl;
   cout << j(500) << endl;
   // Class interpolant knows how to integrate itself exactly.
   cout << i.integral(300 * u::nm, 900 * u::nm) << endl;

   // Product of interpolants is interpolant with no loss of information.
   interpolant<length, double> k("file2.txt", u::nm, 1.0);
   interpolant<length, num::time> m = i * k;

   return 0;
}
```

A linear interpolant can also be initialized from an existing, continuous
function.

```.cpp
#include <cmath> // for erf()
#include <interpolant.hpp>

using namespace num;

int main()
{
   // 'interpolantd' is an alias for 'interpolant<double, double>'.
   interpolantd const i(erf, -1.0, 2.0);

   // Now you can use, for example i(1.3), to do a quick table-lookup of erf(),
   // etc.  This is useful especially for a continuous function that takes a
   // long time to calculate, when it has to be used over and over again in a
   // loop.

   return 0;
}
```

