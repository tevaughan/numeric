
// Copyright 2016  Thomas E. Vaughan
//
// This software is distributable under the terms of the GNU LGPL, Version 3 or
// later.

#include <cmath> // for erf()

#include "catch.hpp"
#include "integral.hpp"
#include "interpolant.hpp"
#include "units.hpp"

using namespace GiNaC;
using namespace num;
using namespace num::u;
using namespace std;

TEST_CASE("Verify interpolation on vector of points.", "[interpolant]")
{
   ilist<double, double> list = {{0.00, 0.00}, {0.50, 0.25}, {1.00, 1.00}};

   auto j = make_linear_interp(list);

   REQUIRE(j(-1.0) == 0.0);
   REQUIRE(j(0.0) == 0.0);
   REQUIRE(dbl(j(0.25)) == Approx(0.25 / 2.0));
   REQUIRE(j(0.50) == 0.25);
   REQUIRE(dbl(j(0.75)) == Approx(1.25 / 2.0));
   REQUIRE(j(1.00) == 1.00);
   REQUIRE(j(1.01) == 0.00); // sparse_table is zero outside bounds.
   REQUIRE(j(2.00) == 0.00); // sparse_table is zero outside bounds.
}

TEST_CASE("Verify integral of interpolant.", "[interpolant]")
{
   ilist<double, double> list = {{0.00, 0.50}, {0.50, 1.00}, {1.00, -1.00}};

   auto i2 = make_linear_interp(list);

   REQUIRE(i2.integral(-0.50, -0.10) == 0.0);
   REQUIRE(i2.integral(-0.25, +0.25) == 0.25 * 0.5 * (0.5 + 0.75));
   REQUIRE(i2.integral(+0.25, +1.25) == 0.5 * 0.25 * (0.75 + 1.0));
   REQUIRE(dbl(i2.integral(+0.55, +0.95)) == Approx(0.0));
   REQUIRE(i2.integral(+0.75, +1.25) == -0.125);
   REQUIRE(i2.integral(+1.25, +0.75) == +0.125);
   REQUIRE(i2.integral(+1.25, +1.50) == 0.00);
   REQUIRE(i2.integral(+1.50, +1.25) == 0.00);
}

TEST_CASE("Verify product of interpolants.", "[interpolant]")
{
   ilist<double, double> list1 = {{0.00, 0.00}, {0.50, 0.25}, {1.00, 1.00}};
   ilist<double, double> list2 = {{0.25, 2.00}, {0.75, 1.00}, {1.25, 0.00}};

   auto i3 = make_linear_interp(list1);
   auto i4 = make_linear_interp(list2);

   auto j3 = i3 * i4;
   auto j4 = i4 * i3;

   REQUIRE(j3(-1.0) == 0.0);
   REQUIRE(j3(0.0) == 0.0);
   REQUIRE(j3(0.25) == 0.25);
   REQUIRE(j3(0.50) == 0.25 * 1.5);
   REQUIRE(dbl(j3(0.75)) == Approx(1.25 / 2.0));
   REQUIRE(j3(1.00) == 0.5);
   REQUIRE(j3(2.00) == 0.0);

   REQUIRE(j4(-1.0) == 0.0);
   REQUIRE(j4(0.0) == 0.0);
   REQUIRE(j4(0.25) == 0.25);
   REQUIRE(j4(0.50) == 0.25 * 1.5);
   REQUIRE(dbl(j4(0.75)) == Approx(1.25 / 2.0));
   REQUIRE(j4(1.00) == 0.5);
   REQUIRE(j4(2.00) == 0.0);
}

TEST_CASE("Verify quotient of interpolants.", "[interpolant]")
{
   ilist<double, double> list1 = {{0.00, 0.00}, {0.50, 0.25}, {1.00, 1.00}};
   ilist<double, double> list2 = {{0.25, 2.00}, {0.75, 1.00}, {1.25, 0.50}};

   auto i3 = make_linear_interp(list1);
   auto i4 = make_linear_interp(list2);

   auto j2 = i3 / i4;

   REQUIRE(j2(-1.0) == 0.0);
   REQUIRE(j2(0.0) == 0.0);
   REQUIRE(j2(0.25) == (0.25 / 2.0) / 2.0);
   REQUIRE(dbl(j2(0.50)) == Approx(0.25 / 1.5));
   REQUIRE(dbl(j2(0.75)) == Approx(1.25 / 2.0));
   REQUIRE(dbl(j2(1.00)) == Approx(1.0 / 0.75));
}

TEST_CASE("Verify file input.", "[interpolant]")
{
   REQUIRE_THROWS(make_linear_interp("nonexistent-file"));
   REQUIRE_THROWS(make_linear_interp("interpolant_test-badinput.txt"));
   REQUIRE_THROWS(make_linear_interp("interpolant_test-badinput2.txt"));

   auto j = make_linear_interp("interpolant_test.txt");

   REQUIRE(j(-1) == 0.0);
   REQUIRE(j(299) == 0.0);
   REQUIRE(dbl(j(300)) == Approx(0.0));
   REQUIRE(dbl(j(350)) == Approx(40.0));
   REQUIRE(j(400) == 80.0);
   REQUIRE(j(450) == 80.0);
   REQUIRE(j(500) == 80.0);
   REQUIRE(dbl(j(900)) == Approx(80.0));
   REQUIRE(dbl(j(950)) == Approx(40.0));
   REQUIRE(dbl(j(1000)) == Approx(0.0));
   REQUIRE(j(1100) == 0.0);

   REQUIRE_THROWS(make_linear_interp("interpolant_test-shortinput.txt"));
}

TEST_CASE("Verify file input and interoperation with units.", "[interpolant]")
{
   using time = num::time;
   auto jj = make_linear_interp<dyndim, dyndim>("interpolant_test.txt", nm, s);
   auto kk = make_linear_interp<length, time>("interpolant_test.txt", nm, s);

   REQUIRE(jj(-1 * nm) == 0.0 * s);
   REQUIRE(jj(299 * nm) == 0.0 * s);
   REQUIRE(dbl(jj(300 * nm) / s) == Approx(0.0));
   REQUIRE(dbl(jj(350 * nm) / s) == Approx(40.0));
   REQUIRE(dbl(jj(400 * nm) / s) == Approx(80.0));
   REQUIRE(jj(450 * nm) == 80.0 * s);
   REQUIRE(jj(500 * nm) == 80.0 * s);
   REQUIRE(jj(900 * nm) == 80.0 * s);
   REQUIRE(dbl(jj(950 * nm) / s) == Approx(40.0));
   REQUIRE(dbl(jj(1000 * nm) / s) == Approx(0.0));
   REQUIRE(jj(1100 * nm) == 0.0 * s);

   REQUIRE(kk(-1 * nm) == 0.0 * s);
   REQUIRE(kk(299 * nm) == 0.0 * s);
   REQUIRE(dbl(kk(300 * nm) / s) == Approx(0.0));
   REQUIRE(dbl(kk(350 * nm) / s) == Approx(40.0));
   REQUIRE(dbl(kk(400 * nm) / s) == Approx(80.0));
   REQUIRE(kk(450 * nm) == 80.0 * s);
   REQUIRE(kk(500 * nm) == 80.0 * s);
   REQUIRE(kk(900 * nm) == 80.0 * s);
   REQUIRE(dbl(kk(950 * nm) / s) == Approx(40.0));
   REQUIRE(dbl(kk(1000 * nm) / s) == Approx(0.0));
   REQUIRE(kk(1100 * nm) == 0.0 * s);
}

TEST_CASE("Verify interoperation with units.", "[interpolant]")
{
   ilist<length, double> list1 = {{300 * nm, 0.0},
                                  {400 * nm, 15.0},
                                  {500 * nm, 80.0},
                                  {900 * nm, 80.0},
                                  {1000 * nm, 0.0}};

   ilist<dyndim, dyndim> list2 = {{300 * nm, 0.0 * s},
                                  {400 * nm, 15.0 * s},
                                  {500 * nm, 80.0 * s},
                                  {900 * nm, 80.0 * s},
                                  {1000 * nm, 0.0 * s}};

   auto i = make_linear_interp<length, double>(list1);
   auto j = make_linear_interp<dyndim, dyndim>(list2);

   REQUIRE(i(-1 * nm) == 0.0);
   REQUIRE(i(299 * nm) == 0.0);
   REQUIRE(dbl(i(300 * nm)) == Approx(0.0));
   REQUIRE(dbl(i(350 * nm)) == Approx(7.5));
   REQUIRE(dbl(i(400 * nm)) == Approx(15.0));
   REQUIRE(dbl(i(450 * nm)) == Approx(0.5 * (15.0 + 80.0)));
   REQUIRE(i(500 * nm) == 80.0);
   REQUIRE(i(600 * nm) == 80.0);
   REQUIRE(dbl(i(900 * nm)) == Approx(80.0));
   REQUIRE(dbl(i(950 * nm)) == Approx(40.0));
   REQUIRE(dbl(i(1000 * nm)) == Approx(0.0));
   REQUIRE(i(1100 * nm) == 0.0);

   REQUIRE(j(-1 * nm) == 0.0 * s);
   REQUIRE(j(299 * nm) == 0.0 * s);
   REQUIRE(dbl(j(300 * nm) / s) == Approx(0.0));
   REQUIRE(dbl(j(350 * nm) / s) == Approx(7.5));
   REQUIRE(dbl(j(400 * nm) / s) == Approx(15.0));
   REQUIRE(dbl(j(450 * nm) / s) == Approx(0.5 * (15.0 + 80.0)));
   REQUIRE(j(500 * nm) == 80.0 * s);
   REQUIRE(j(600 * nm) == 80.0 * s);
   REQUIRE(dbl(j(900 * nm) / s) == Approx(80.0));
   REQUIRE(dbl(j(950 * nm) / s) == Approx(40.0));
   REQUIRE(dbl(j(1000 * nm) / s) == Approx(0.0));
   REQUIRE(j(1100 * nm) == 0.0 * s);
}

TEST_CASE(
      "Verify scalar multiplication of interpolant on right.", "[interpolant]")
{
   using time = num::time;
   auto i3 = make_linear_interp<length, time>("interpolant_test.txt", nm, s);
   auto i4 = make_linear_interp<dyndim, dyndim>("interpolant_test.txt", nm, s);
   auto j3 = i3 * (2 * m / s);
   auto j4 = i4 * (2 * m / s);

   REQUIRE(j3(-1 * nm) == 0.0 * m);
   REQUIRE(j3(299 * nm) == 0.0 * m);
   REQUIRE(dbl(j3(300 * nm) / m) == Approx(0.0));
   REQUIRE(dbl(j3(350 * nm) / m) == Approx(80.0));
   REQUIRE(dbl(j3(400 * nm) / m) == Approx(160.0));
   REQUIRE(j3(450 * nm) == 160.0 * m);
   REQUIRE(j3(500 * nm) == 160.0 * m);
   REQUIRE(j3(900 * nm) == 160.0 * m);
   REQUIRE(dbl(j3(950 * nm) / m) == Approx(80.0));
   REQUIRE(dbl(j3(1000 * nm) / m) == Approx(0.0));
   REQUIRE(j3(1100 * nm) == 0.0 * m);

   REQUIRE(j4(-1 * nm) == 0.0 * m);
   REQUIRE(j4(299 * nm) == 0.0 * m);
   REQUIRE(dbl(j4(300 * nm) / m) == Approx(0.0));
   REQUIRE(dbl(j4(350 * nm) / m) == Approx(80.0));
   REQUIRE(dbl(j4(400 * nm) / m) == Approx(160.0));
   REQUIRE(j4(450 * nm) == 160.0 * m);
   REQUIRE(j4(500 * nm) == 160.0 * m);
   REQUIRE(j4(900 * nm) == 160.0 * m);
   REQUIRE(dbl(j4(950 * nm) / m) == Approx(80.0));
   REQUIRE(dbl(j4(1000 * nm) / m) == Approx(0.0));
   REQUIRE(j4(1100 * nm) == 0.0 * m);

   auto k3 = i3 * 2.0;
   auto k4 = i4 * 2.0;

   REQUIRE(k3(-1 * nm) == 0.0 * s);
   REQUIRE(k3(299 * nm) == 0.0 * s);
   REQUIRE(dbl(k3(300 * nm) / s) == Approx(0.0));
   REQUIRE(dbl(k3(350 * nm) / s) == Approx(80.0));
   REQUIRE(dbl(k3(400 * nm) / s) == Approx(160.0));
   REQUIRE(k3(450 * nm) == 160.0 * s);
   REQUIRE(k3(500 * nm) == 160.0 * s);
   REQUIRE(k3(900 * nm) == 160.0 * s);
   REQUIRE(dbl(k3(950 * nm) / s) == Approx(80.0));
   REQUIRE(dbl(k3(1000 * nm) / s) == Approx(0.0));
   REQUIRE(k3(1100 * nm) == 0.0 * s);

   REQUIRE(k4(-1 * nm) == 0.0 * s);
   REQUIRE(k4(299 * nm) == 0.0 * s);
   REQUIRE(dbl(k4(300 * nm) / s) == Approx(0.0));
   REQUIRE(dbl(k4(350 * nm) / s) == Approx(80.0));
   REQUIRE(dbl(k4(400 * nm) / s) == Approx(160.0));
   REQUIRE(k4(450 * nm) == 160.0 * s);
   REQUIRE(k4(500 * nm) == 160.0 * s);
   REQUIRE(k4(900 * nm) == 160.0 * s);
   REQUIRE(dbl(k4(950 * nm) / s) == Approx(80.0));
   REQUIRE(dbl(k4(1000 * nm) / s) == Approx(0.0));
   REQUIRE(k4(1100 * nm) == 0.0 * s);

   i3 *= 2.0;
   i4 *= 2.0;

   REQUIRE(i3(-1 * nm) == 0.0 * s);
   REQUIRE(i3(299 * nm) == 0.0 * s);
   REQUIRE(dbl(i3(300 * nm) / s) == Approx(0.0));
   REQUIRE(dbl(i3(350 * nm) / s) == Approx(80.0));
   REQUIRE(dbl(i3(400 * nm) / s) == Approx(160.0));
   REQUIRE(i3(450 * nm) == 160.0 * s);
   REQUIRE(i3(500 * nm) == 160.0 * s);
   REQUIRE(i3(900 * nm) == 160.0 * s);
   REQUIRE(dbl(i3(950 * nm) / s) == Approx(80.0));
   REQUIRE(dbl(i3(1000 * nm) / s) == Approx(0.0));
   REQUIRE(i3(1100 * nm) == 0.0 * s);

   REQUIRE(i4(-1 * nm) == 0.0 * s);
   REQUIRE(i4(299 * nm) == 0.0 * s);
   REQUIRE(dbl(i4(300 * nm) / s) == Approx(0.0));
   REQUIRE(dbl(i4(350 * nm) / s) == Approx(80.0));
   REQUIRE(dbl(i4(400 * nm) / s) == Approx(160.0));
   REQUIRE(i4(450 * nm) == 160.0 * s);
   REQUIRE(i4(500 * nm) == 160.0 * s);
   REQUIRE(i4(900 * nm) == 160.0 * s);
   REQUIRE(dbl(i4(950 * nm) / s) == Approx(80.0));
   REQUIRE(dbl(i4(1000 * nm) / s) == Approx(0.0));
   REQUIRE(i4(1100 * nm) == 0.0 * s);
}

TEST_CASE("Verify scalar multiplication on left.", "[interpolant]")
{
   using time = num::time;
   auto i3 = make_linear_interp<length, time>("interpolant_test.txt", nm, s);
   auto i4 = make_linear_interp<dyndim, dyndim>("interpolant_test.txt", nm, s);
   auto j3 = (2 * m / s) * i3;
   auto j4 = (2 * m / s) * i4;

   REQUIRE(j3(-1 * nm) == 0.0 * m);
   REQUIRE(j3(299 * nm) == 0.0 * m);
   REQUIRE(dbl(j3(300 * nm) / m) == Approx(0.0));
   REQUIRE(dbl(j3(350 * nm) / m) == Approx(80.0));
   REQUIRE(dbl(j3(400 * nm) / m) == Approx(160.0));
   REQUIRE(j3(450 * nm) == 160.0 * m);
   REQUIRE(j3(500 * nm) == 160.0 * m);
   REQUIRE(j3(900 * nm) == 160.0 * m);
   REQUIRE(dbl(j3(950 * nm) / m) == Approx(80.0));
   REQUIRE(dbl(j3(1000 * nm) / m) == Approx(0.0));
   REQUIRE(j3(1100 * nm) == 0.0 * m);

   REQUIRE(j4(-1 * nm) == 0.0 * m);
   REQUIRE(j4(299 * nm) == 0.0 * m);
   REQUIRE(dbl(j4(300 * nm) / m) == Approx(0.0));
   REQUIRE(dbl(j4(350 * nm) / m) == Approx(80.0));
   REQUIRE(dbl(j4(400 * nm) / m) == Approx(160.0));
   REQUIRE(j4(450 * nm) == 160.0 * m);
   REQUIRE(j4(500 * nm) == 160.0 * m);
   REQUIRE(j4(900 * nm) == 160.0 * m);
   REQUIRE(dbl(j4(950 * nm) / m) == Approx(80.0));
   REQUIRE(dbl(j4(1000 * nm) / m) == Approx(0.0));
   REQUIRE(j4(1100 * nm) == 0.0 * m);

   auto k3 = 2.0 * i3;
   auto k4 = 2.0 * i4;

   REQUIRE(k3(-1 * nm) == 0.0 * s);
   REQUIRE(k3(299 * nm) == 0.0 * s);
   REQUIRE(dbl(k3(300 * nm) / s) == Approx(0.0));
   REQUIRE(dbl(k3(350 * nm) / s) == Approx(80.0));
   REQUIRE(dbl(k3(400 * nm) / s) == Approx(160.0));
   REQUIRE(k3(450 * nm) == 160.0 * s);
   REQUIRE(k3(500 * nm) == 160.0 * s);
   REQUIRE(k3(900 * nm) == 160.0 * s);
   REQUIRE(dbl(k3(950 * nm) / s) == Approx(80.0));
   REQUIRE(dbl(k3(1000 * nm) / s) == Approx(0.0));
   REQUIRE(k3(1100 * nm) == 0.0 * s);

   REQUIRE(k4(-1 * nm) == 0.0 * s);
   REQUIRE(k4(299 * nm) == 0.0 * s);
   REQUIRE(dbl(k4(300 * nm) / s) == Approx(0.0));
   REQUIRE(dbl(k4(350 * nm) / s) == Approx(80.0));
   REQUIRE(dbl(k4(400 * nm) / s) == Approx(160.0));
   REQUIRE(k4(450 * nm) == 160.0 * s);
   REQUIRE(k4(500 * nm) == 160.0 * s);
   REQUIRE(k4(900 * nm) == 160.0 * s);
   REQUIRE(dbl(k4(950 * nm) / s) == Approx(80.0));
   REQUIRE(dbl(k4(1000 * nm) / s) == Approx(0.0));
   REQUIRE(k4(1100 * nm) == 0.0 * s);
}

TEST_CASE("Verify scalar division of interpolant.", "[interpolant]")
{
   using time = num::time;
   auto i3 = make_linear_interp<length, time>("interpolant_test.txt", nm, s);
   auto i4 = make_linear_interp<dyndim, dyndim>("interpolant_test.txt", nm, s);
   auto j3 = i3 / (2 * s);
   auto j4 = i4 / (2 * s);

   REQUIRE(j3(-1 * nm) == 0.0);
   REQUIRE(j3(299 * nm) == 0.0);
   REQUIRE(dbl(j3(300 * nm)) == Approx(0.0));
   REQUIRE(dbl(j3(350 * nm)) == Approx(20.0));
   REQUIRE(dbl(j3(400 * nm)) == Approx(40.0));
   REQUIRE(j3(450 * nm) == 40.0);
   REQUIRE(j3(500 * nm) == 40.0);
   REQUIRE(j3(900 * nm) == 40.0);
   REQUIRE(dbl(j3(950 * nm)) == Approx(20.0));
   REQUIRE(dbl(j3(1000 * nm)) == Approx(0.0));
   REQUIRE(j3(1100 * nm) == 0.0);

   REQUIRE(dbl(j4(-1 * nm)) == 0.0);
   REQUIRE(dbl(j4(299 * nm)) == 0.0);
   REQUIRE(dbl(j4(300 * nm)) == Approx(0.0));
   REQUIRE(dbl(j4(350 * nm)) == Approx(20.0));
   REQUIRE(dbl(j4(400 * nm)) == Approx(40.0));
   REQUIRE(dbl(j4(450 * nm)) == 40.0);
   REQUIRE(dbl(j4(500 * nm)) == 40.0);
   REQUIRE(dbl(j4(900 * nm)) == 40.0);
   REQUIRE(dbl(j4(950 * nm)) == Approx(20.0));
   REQUIRE(dbl(j4(1000 * nm)) == Approx(0.0));
   REQUIRE(dbl(j4(1100 * nm)) == 0.0);

   i3 /= 2.0;
   i4 /= 2.0;

   REQUIRE(i3(-1 * nm) == 0.0 * s);
   REQUIRE(i3(299 * nm) == 0.0 * s);
   REQUIRE(dbl(i3(300 * nm) / s) == Approx(0.0));
   REQUIRE(dbl(i3(350 * nm) / s) == Approx(20.0));
   REQUIRE(dbl(i3(400 * nm) / s) == Approx(40.0));
   REQUIRE(i3(450 * nm) == 40.0 * s);
   REQUIRE(i3(500 * nm) == 40.0 * s);
   REQUIRE(i3(900 * nm) == 40.0 * s);
   REQUIRE(dbl(i3(950 * nm) / s) == Approx(20.0));
   REQUIRE(dbl(i3(1000 * nm) / s) == Approx(0.0));
   REQUIRE(i3(1100 * nm) == 0.0 * s);

   REQUIRE(i4(-1 * nm) == 0.0 * s);
   REQUIRE(i4(299 * nm) == 0.0 * s);
   REQUIRE(dbl(i4(300 * nm) / s) == Approx(0.0));
   REQUIRE(dbl(i4(350 * nm) / s) == Approx(20.0));
   REQUIRE(dbl(i4(400 * nm) / s) == Approx(40.0));
   REQUIRE(i4(450 * nm) == 40.0 * s);
   REQUIRE(i4(500 * nm) == 40.0 * s);
   REQUIRE(i4(900 * nm) == 40.0 * s);
   REQUIRE(dbl(i4(950 * nm) / s) == Approx(20.0));
   REQUIRE(dbl(i4(1000 * nm) / s) == Approx(0.0));
   REQUIRE(i4(1100 * nm) == 0.0 * s);
}

double my_erf(double x) { return erf(x); }

TEST_CASE(
      "Verify that interpolant of function produces right integral.",
      "[interpolant]")
{
   std::function<area(length)> square = [](length x) { return x * x; };
   volume const i                     = num::integral(square, 0 * cm, 1 * cm);

   auto j3 = make_linear_interp<length, area>(square, 0 * cm, 1 * cm);
   auto j4 = make_linear_interp<dyndim, dyndim>(square, 0 * cm, 1 * cm);

   REQUIRE(dbl(j3.integral() / i) == Approx(1.0));
   REQUIRE(dbl(j3.integral(0 * cm, 1 * cm) / i) == Approx(1.0));

   REQUIRE(dbl(j4.integral() / i) == Approx(1.0));
   REQUIRE(dbl(j4.integral(0 * cm, 1 * cm) / i) == Approx(1.0));

   auto k3 = make_linear_interp<length, area>(square, 1 * cm, 0 * cm);
   auto k4 = make_linear_interp<dyndim, dyndim>(square, 1 * cm, 0 * cm);

   REQUIRE(dbl(k3.integral(0 * cm, 1 * cm) / i) == Approx(1.0));
   REQUIRE(dbl(k4.integral(0 * cm, 1 * cm) / i) == Approx(1.0));

   REQUIRE_THROWS(make_linear_interp(square, 0 * cm, 1 * cm, -1.0E-06));

   double     tol = 1.0E-04;
   auto const e1  = make_linear_interp(erf, -1.0, 2.0, tol);

   REQUIRE(
         dbl(e1.integral(-1.0, 2.0)) / num::integral(my_erf, -1.0, 2.0, tol) ==
         Approx(1.0).epsilon(tol));

   tol                             = 1.0E-03;
   std::function<double(double)> g = [](double x) {
      return exp(-0.5 * x * x);
   };
   auto const ig = make_linear_interp(g, -5.0, +5.0, tol);
   REQUIRE(dbl(ig.integral()) == Approx(sqrt(2.0 * M_PI)).epsilon(tol));
}

