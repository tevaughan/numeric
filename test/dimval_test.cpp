
#include "catch.hpp"
#include "dimval.hpp"

using namespace num;
using namespace std;

TEST_CASE("Verify basic operation of dimval.", "[dimval]")
{
   length const len1 = meters(3.0);
   length const len2 = 6.0 * meter;
   REQUIRE(len1 / len2 == Approx(0.5));
}

