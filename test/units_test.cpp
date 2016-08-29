
#include "catch.hpp"
#include "units.hpp"

using namespace num;
using namespace std;

TEST_CASE("Verify basic operation of units.", "[units]")
{
   length const len1 = meters(3.0);
   length const len2 = 6.0 * m;
   REQUIRE(len1 / len2 == Approx(0.5));
}

