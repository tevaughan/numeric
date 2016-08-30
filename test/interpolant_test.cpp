
#include "catch.hpp"
#include "interpolant.hpp"
#include "units.hpp"

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

TEST_CASE("Verify file input.", "[interpolant]")
{
   interpolantd i("interpolant_test.txt");
   REQUIRE(i(-1) == 0.0);
   REQUIRE(i(299) == 0.0);
   REQUIRE(i(300) == 0.0);
   REQUIRE(i(350) == Approx(40.0));
   REQUIRE(i(400) == 80.0);
   REQUIRE(i(450) == 80.0);
   REQUIRE(i(500) == 80.0);
   REQUIRE(i(900) == 80.0);
   REQUIRE(i(950) == Approx(40.0));
   REQUIRE(i(1000) == 0.0);
   REQUIRE(i(1100) == 0.0);
}

TEST_CASE("Verify file input and interoperation with units.", "[interpolant]")
{
   interpolant<length, num::time> i("interpolant_test.txt", nm, s);
   REQUIRE(i(-1 * nm) == 0.0 * s);
   REQUIRE(i(299 * nm) == 0.0 * s);
   REQUIRE(i(300 * nm) == 0.0 * s);
   REQUIRE(i(350 * nm) / s == Approx(40.0));
   REQUIRE(i(400 * nm) == 80.0 * s);
   REQUIRE(i(450 * nm) == 80.0 * s);
   REQUIRE(i(500 * nm) == 80.0 * s);
   REQUIRE(i(900 * nm) == 80.0 * s);
   REQUIRE(i(950 * nm) / s == Approx(40.0));
   REQUIRE(i(1000 * nm) == 0.0 * s);
   REQUIRE(i(1100 * nm) == 0.0 * s);
}

TEST_CASE("Verify interoperation with units.", "[interpolant]")
{
   ilist<length, double> list = {{300 * nm, 0.0},
                                 {400 * nm, 15.0},
                                 {500 * nm, 80.0},
                                 {900 * nm, 80.0},
                                 {1000 * nm, 0.0}};
   interpolant<length, double> i(list);
   REQUIRE(i(-1 * nm) == 0.0);
   REQUIRE(i(299 * nm) == 0.0);
   REQUIRE(i(300 * nm) == 0.0);
   REQUIRE(i(350 * nm) == Approx(7.5));
   REQUIRE(i(400 * nm) == 15.0);
   REQUIRE(i(450 * nm) == Approx(0.5 * (15.0 + 80.0)));
   REQUIRE(i(500 * nm) == 80.0);
   REQUIRE(i(600 * nm) == 80.0);
   REQUIRE(i(900 * nm) == 80.0);
   REQUIRE(i(950 * nm) == Approx(40.0));
   REQUIRE(i(1000 * nm) == 0.0);
   REQUIRE(i(1100 * nm) == 0.0);
}

