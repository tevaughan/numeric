
// Copyright 2016  Thomas E. Vaughan
//
// This software is distributable under the terms of the GNU GPL, Version 3 or
// later.

#include <sstream> // for ostringstream

#include "catch.hpp"
#include "integral.hpp"
#include "interpolant.hpp"
#include "units.hpp"

using namespace num;
using namespace std;

TEST_CASE("Verify integration of lambda.", "[integral]")
{
   function<double(double)> linear = [](double x) { return x; };
   // Verify that integral of x of 0 to 1 is 1/2.
   REQUIRE(integral(linear, 0, 1) == Approx(0.5));
}

area square(length x) { return x * x; }

TEST_CASE("Verify integration of (dimval) global function.", "[integral]")
{
   volume const i = integral(square, 0 * cm, 1 * cm);
   ostringstream oss;
   oss << i;
   // Verify that integral of x^2 from  0 cm  to  1 cm  is  1/3 cm^3.
   REQUIRE(oss.str() == "[3.33333e-07 m^3]");
   REQUIRE(i / cm.pow<3>() == Approx(1.0 / 3.0));
}

struct my_sin {
   double freq;
   double sin(double x) { return std::sin(freq * x); }
   explicit my_sin(double f = 1.0) : freq(f) {}
};

TEST_CASE("Verify integration of member function (via lambda).", "[integral]")
{
   my_sin f;
   function<double(double)> s = [&f](double x) { return f.sin(x); };
   double const i = integral(s, 0, M_PI);
   // Verify that integral of sin(x) from 0 to pi is 2.
   REQUIRE(i == Approx(2.0));
   f.freq = 2.0;
   double const j = integral(s, 0, M_PI);
   // Verify that integral of sin(2*x) from 0 to pi is 0.
   REQUIRE(j == Approx(0.0));
}

TEST_CASE("Verify limit of tolerance.", "[integral]")
{
   my_sin f;
   function<double(double)> s = [&f](double x) { return f.sin(x); };
   double constexpr tol = 1.0E-11;
   double const i = integral(s, 0, M_PI, tol);
   // Verify that integral of sin(x) from 0 to pi is 2.
   REQUIRE(i == Approx(2.0));
   f.freq = 2.0;
   double const j = integral(s, 0, M_PI, tol);
   // Verify that integral of sin(2*x) from 0 to pi is 0.
   REQUIRE(j == Approx(0.0));
}

TEST_CASE("Trigger coverage of code requiring at least two samples.",
          "[integral]")
{
   volume const i = integral(square, 1 * cm, 2 * cm, 1.0E-06, 0);
   // Verify that integral of x^2 from 1 cm  to  2 cm  is  7/3 cm^3.
   REQUIRE(i / cm.pow<3>() == Approx(7.0 / 3.0));
}

// This is a bit silly because interpolant has its own, optimized
// interpolant::integral() member function, but global integral() is tested
// here for completeness.
TEST_CASE("Verify integration of interpolant.", "[integral]")
{
   ilist<double, double> list = {{0.00, 0.50}, {0.50, 1.00}, {1.00, -1.00}};
   interpolantd i(list);
   function<double(double)> f(i);
   REQUIRE(integral(f, -0.50, -0.10) == Approx(0.2));
   REQUIRE(integral(f, -0.25, +0.25) ==
           Approx(0.25 * 0.5 + 0.25 * 0.5 * (0.5 + 0.75)));
   REQUIRE(integral(f, +0.25, +1.25) ==
           Approx(0.5 * 0.25 * (0.75 + 1.0) - 0.25));
   REQUIRE(integral(f, +0.55, +0.95) == Approx(0.0));
   REQUIRE(integral(f, +0.75, +1.25) == Approx(-0.375));
   REQUIRE(integral(f, +1.25, +1.50) == Approx(-0.25));
   REQUIRE(integral(f, +1.50, +1.25) == Approx(+0.25));
}

TEST_CASE("Verify throw on illegal tolerance.", "[integral]")
{
   REQUIRE_THROWS(integral(square, 1 * cm, 2 * cm, -1.0E-06));
}

