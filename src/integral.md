
# Integration {#integrals}



```cpp
#include <iostream>
#include <integral.hpp>
#include <units.hpp>

using namespace num;
using namespace std;

area square(length x) { return x * x; }

int main()
{
   volume const i = integral(square, 0 * u::cm, 1 * u::cm);
   cout << i << endl; // "[3.33333e-07 m^3]"
   return 0;
}
```

