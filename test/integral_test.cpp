
#include <sstream> // for ostringstream

#include "catch.hpp"
#include "integral.hpp"
#include "units.hpp"

using namespace num;
using namespace std;

TEST_CASE("Verify numeric integration of linear function.", "[integral]")
{
   auto linear = [](double x) { return x; };
   REQUIRE(integral(linear, 0.0, 1.0) == Approx(0.5));
}

TEST_CASE("Verify interoperation with dimval.", "[integral]")
{
   auto square = [](length x) { return x * x; };
   volume const i = integral(square, 0.0 * cm, 1.0 * cm);
   ostringstream oss;
   oss << i;
   REQUIRE(oss.str() == "[3.33333e-07 m^3]");
   double const j = i / cm.pow<3>();
   double constexpr third = 1.0 / 3.0;
   REQUIRE(j == Approx(third));
}

