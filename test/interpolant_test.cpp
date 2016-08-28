
#include "catch.hpp"
#include "interpolant.hpp"

using namespace num;
using namespace std;

TEST_CASE("Verify interpolation on vector of points.", "[interpolant]")
{
   ilist<double, double> list = {{0.00, 0.00}, {0.50, 0.25}, {1.00, 1.00}};
   interpolantd i(list);
   REQUIRE(i(-1.0) == 0.0);
   REQUIRE(i(0.0) == 0.0);
   REQUIRE(i(0.25) == Approx(0.25/2.0));
   REQUIRE(i(0.50) == 0.25);
   REQUIRE(i(0.75) == Approx(1.25/2.0));
   REQUIRE(i(1.00) == 1.00);
   REQUIRE(i(2.00) == 1.00);
}

