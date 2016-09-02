
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
   REQUIRE(i(0.25) == Approx(0.25 / 2.0));
   REQUIRE(i(0.50) == 0.25);
   REQUIRE(i(0.75) == Approx(1.25 / 2.0));
   REQUIRE(i(1.00) == 1.00);
   REQUIRE(i(2.00) == 1.00);
}

TEST_CASE("Verify product of interpolants.", "[interpolant]")
{
   ilist<double, double> list1 = {{0.00, 0.00}, {0.50, 0.25}, {1.00, 1.00}};
   ilist<double, double> list2 = {{0.25, 2.00}, {0.75, 1.00}, {1.25, 0.00}};
   interpolantd i1(list1);
   interpolantd i2(list2);
   interpolantd i3(i1 * i2);
   REQUIRE(i3(-1.0) == 0.0);
   REQUIRE(i3(0.0) == 0.0);
   REQUIRE(i3(0.25) == 0.25);
   REQUIRE(i3(0.50) == 0.25 * 1.5);
   REQUIRE(i3(0.75) == Approx(1.25 / 2.0));
   REQUIRE(i3(1.00) == 0.5);
   REQUIRE(i3(2.00) == 0.0);
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

TEST_CASE("Verify scalar multiplication of interpolant on right.",
          "[interpolant]")
{
   interpolant<length, num::time> i("interpolant_test.txt", nm, s);
   interpolant<length, length> j = i * (2 * m / s);
   REQUIRE(j(-1 * nm) == 0.0 * m);
   REQUIRE(j(299 * nm) == 0.0 * m);
   REQUIRE(j(300 * nm) == 0.0 * m);
   REQUIRE(j(350 * nm) / m == Approx(80.0));
   REQUIRE(j(400 * nm) == 160.0 * m);
   REQUIRE(j(450 * nm) == 160.0 * m);
   REQUIRE(j(500 * nm) == 160.0 * m);
   REQUIRE(j(900 * nm) == 160.0 * m);
   REQUIRE(j(950 * nm) / m == Approx(80.0));
   REQUIRE(j(1000 * nm) == 0.0 * m);
   REQUIRE(j(1100 * nm) == 0.0 * m);
}

TEST_CASE("Verify scalar multiplication of interpolant on left.",
          "[interpolant]")
{
   interpolant<length, num::time> i("interpolant_test.txt", nm, s);
   interpolant<length, length> j = (2 * m / s) * i;
   REQUIRE(j(-1 * nm) == 0.0 * m);
   REQUIRE(j(299 * nm) == 0.0 * m);
   REQUIRE(j(300 * nm) == 0.0 * m);
   REQUIRE(j(350 * nm) / m == Approx(80.0));
   REQUIRE(j(400 * nm) == 160.0 * m);
   REQUIRE(j(450 * nm) == 160.0 * m);
   REQUIRE(j(500 * nm) == 160.0 * m);
   REQUIRE(j(900 * nm) == 160.0 * m);
   REQUIRE(j(950 * nm) / m == Approx(80.0));
   REQUIRE(j(1000 * nm) == 0.0 * m);
   REQUIRE(j(1100 * nm) == 0.0 * m);
}

TEST_CASE("Verify scalar division of interpolant.", "[interpolant]")
{
   interpolant<length, num::time> i("interpolant_test.txt", nm, s);
   interpolant<length, double> j = i / (2 * s);
   REQUIRE(j(-1 * nm) == 0.0);
   REQUIRE(j(299 * nm) == 0.0);
   REQUIRE(j(300 * nm) == 0.0);
   REQUIRE(j(350 * nm) == Approx(20.0));
   REQUIRE(j(400 * nm) == 40.0);
   REQUIRE(j(450 * nm) == 40.0);
   REQUIRE(j(500 * nm) == 40.0);
   REQUIRE(j(900 * nm) == 40.0);
   REQUIRE(j(950 * nm) == Approx(20.0));
   REQUIRE(j(1000 * nm) == 0.0);
   REQUIRE(j(1100 * nm) == 0.0);
}

