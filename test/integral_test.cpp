
// Copyright 2016-2017  Thomas E. Vaughan
//
// This software is distributable under the terms of the GNU LGPL, Version 3 or
// later.

#include <sstream> // for ostringstream

#include "catch.hpp"
#include "integral.hpp"
#include "interpolant.hpp"
#include "rk.hpp"
#include "units.hpp"

using namespace num;
using namespace num::u;
using namespace std;

TEST_CASE("Verify integration of lambda.", "[integral]")
{
   function<double(double)> linear = [](double x) { return x; };
   // Verify that integral of x from 0 to 1 is 1/2.
   REQUIRE(integral(linear, 0.0, 1.0) == Approx(0.5));
   REQUIRE(rk_quadd(linear, 0.0, 1.0).def_int() == Approx(0.5));
   function<double(double)> f = [](double x) { return 1.0 / (1.0 + x * x); };
   double const max = sqrt(1.0 / numeric_limits<double>::min());
   REQUIRE(integral(f, -max, +max) == Approx(M_PI));
   REQUIRE(rk_quadd(f, -max, +max).def_int() == Approx(M_PI));
}

area square1(length x) { return x * x; }
area square2(dyndim x) { return x * x; }
dyndim square3(length x) { return x * x; }
dyndim square4(dyndim const &x) { return x * x; }

TEST_CASE("Verify integration of (dimval) global function.", "[integral]")
{
   volume const i1 = integral(square1, 0 * cm, 1 * cm);
   volume const i2 = integral(square2, 0 * cm, 1 * cm);
   volume const j1 =
         rk_quad<length, volume>(square1, 0 * cm, 1 * cm).def_int();
   volume const j2 =
         rk_quad<length, volume>(square2, 0 * cm, 1 * cm).def_int();
   ostringstream oss1, oss2;
   oss1 << i1;
   oss2 << i2;
   // Verify that integral of x^2 from  0 cm  to  1 cm  is  1/3 cm^3.
   REQUIRE(oss1.str() == "[3.33333e-07 m^3]");
   REQUIRE(oss2.str() == "[3.33333e-07 m^3]");
   REQUIRE(i1 / pow<3>(cm) == Approx(1.0 / 3.0));
   REQUIRE(i2 / pow<3>(cm) == Approx(1.0 / 3.0));
   REQUIRE(j1 / pow<3>(cm) == Approx(1.0 / 3.0));
   REQUIRE(j2 / pow<3>(cm) == Approx(1.0 / 3.0));
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
   double const i = integral(s, 0, M_PI, 1.0E-09);
   double const k = rk_quadd(s, 0, M_PI, 1.0E-09).def_int();
   // Verify that integral of sin(x) from 0 to pi is 2.
   REQUIRE(i == Approx(2.0));
   REQUIRE(k == Approx(2.0));
   f.freq = 2.0;
   double const j = integral(s, 0, M_PI);
   double const m = rk_quadd(s, 0, M_PI).def_int();
   // Verify that integral of sin(2*x) from 0 to pi is 0.
   REQUIRE(j == Approx(0.0));
   REQUIRE(m == Approx(0.0));
}

double catch_epsilon(double tol)
{
   double constexpr eps = 100.0 * numeric_limits<double>::epsilon();
   return (tol < eps ? eps : tol);
}

TEST_CASE("Verify limit of tolerance.", "[integral]")
{
   my_sin f;
   function<double(double)> s = [&f](double x) { return f.sin(x); };
   double constexpr tol = 1.0E-17;
   double const i = integral(s, 0, M_PI, tol);
   double const k = rk_quadd(s, 0, M_PI, tol).def_int();
   // Verify that integral of sin(x) from 0 to pi is 2.
   REQUIRE(i == Approx(2.0).epsilon(10000 * tol));
   REQUIRE(k == Approx(2.0).epsilon(10000 * tol));
   f.freq = 2.0;
   double const j = integral(s, 0, M_PI, tol);
   double const m = rk_quadd(s, 0, M_PI, tol).def_int();
   // Verify that integral of sin(2*x) from 0 to pi is 0.
   REQUIRE(j == Approx(0.0).epsilon(10000 * tol));
   REQUIRE(m == Approx(0.0).epsilon(10000 * tol));
}

TEST_CASE("Trigger coverage of code requiring at least two samples.",
          "[integral]")
{
   volume const i = integral(square1, 1 * cm, 2 * cm, 1.0E-06, 0);
   volume const j = rk_quad<length, volume>(square1, 1 * cm, 2 * cm, 1.0E-06,
                                            0).def_int();
   // Verify that integral of x^2 from 1 cm  to  2 cm  is  7/3 cm^3.
   REQUIRE(i / pow<3>(cm) == Approx(7.0 / 3.0));
   REQUIRE(j / pow<3>(cm) == Approx(7.0 / 3.0));
}

// This is a bit silly because interpolant has its own, optimized
// interpolant::integral() member function, but global integral() is tested
// here for completeness.
TEST_CASE("Verify integration of interpolant.", "[integral]")
{
   ilist<double, double> list = {{0.00, 0.50}, {0.50, 1.00}, {1.00, -1.00}};
   auto i = make_linear_interp(list);
   function<double(double)> f = [&i](double x) {
      return GiNaC::ex_to<GiNaC::numeric>(i(x)).to_double();
   };
   REQUIRE(integral(f, -0.50, -0.10) == Approx(0.0));
   REQUIRE(rk_quadd(f, -0.50, -0.10).def_int() == Approx(0.0));
   REQUIRE(integral(f, -0.25, +0.25) ==
           Approx(0.25 * 0.5 * (0.5 + 0.75)));
   REQUIRE(rk_quadd(f, -0.25, +0.25).def_int() ==
           Approx(0.25 * 0.5 * (0.5 + 0.75)));
   REQUIRE(integral(f, +0.25, +1.25) ==
           Approx(0.5 * 0.25 * (0.75 + 1.0)));
   REQUIRE(rk_quadd(f, +0.25, +1.25).def_int() ==
           Approx(0.5 * 0.25 * (0.75 + 1.0)));
   REQUIRE(integral(f, +0.55, +0.95) == Approx(0.0));
   REQUIRE(rk_quadd(f, +0.55, +0.95).def_int() == Approx(0.0));
   REQUIRE(integral(f, +0.75, +1.25) == Approx(-0.125));
   REQUIRE(rk_quadd(f, +0.75, +1.25).def_int() == Approx(-0.125));
   REQUIRE(integral(f, +1.25, +1.50) == Approx(0.0));
   REQUIRE(rk_quadd(f, +1.25, +1.50).def_int() == Approx(0.0));
   REQUIRE(integral(f, +1.50, +1.25) == Approx(0.0));
   REQUIRE(rk_quadd(f, +1.50, +1.25).def_int() == Approx(0.0));
}

TEST_CASE("Verify throw on illegal tolerance.", "[integral]")
{
   REQUIRE_THROWS(integral(square1, 1 * cm, 2 * cm, -1.0E-06));
   REQUIRE_THROWS((rk_quad<length, volume>(square1, 1 * cm, 2 * cm, -1.0E-06)
                         .def_int()));
}

