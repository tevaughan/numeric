
// Copyright 2016  Thomas E. Vaughan
//
// This software is distributable under the terms of the GNU LGPL, Version 3 or
// later.

#include <cmath> // for erf()

#include "catch.hpp"
#include "integral.hpp"
#include "interpolant.hpp"
#include "units.hpp"

using namespace num;
using namespace num::u;
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

TEST_CASE("Verify integral of interpolant.", "[interpolant]")
{
   ilist<double, double> list = {{0.00, 0.50}, {0.50, 1.00}, {1.00, -1.00}};
   interpolantd i(list);
   REQUIRE(i.integral(-0.50, -0.10) == 0.2);
   REQUIRE(i.integral(-0.25, +0.25) == 0.25 * 0.5 + 0.25 * 0.5 * (0.5 + 0.75));
   REQUIRE(i.integral(+0.25, +1.25) == 0.5 * 0.25 * (0.75 + 1.0) - 0.25);
   REQUIRE(i.integral(+0.55, +0.95) == 0.0);
   REQUIRE(i.integral(+0.75, +1.25) == -0.375);
   REQUIRE(i.integral(+1.25, +1.50) == -0.25);
   REQUIRE(i.integral(+1.50, +1.25) == +0.25);
   interpolantd j;
   REQUIRE_THROWS(j.integral(0.0, 1.0));
}

TEST_CASE("Verify product of interpolants.", "[interpolant]")
{
   ilist<double, double> list1 = {{0.00, 0.00}, {0.50, 0.25}, {1.00, 1.00}};
   ilist<double, double> list2 = {{0.25, 2.00}, {0.75, 1.00}, {1.25, 0.00}};
   interpolantd i1(list1);
   interpolantd i2(list2);
   interpolantd i3(i1 * i2);
   interpolantd i4(i2 * i1);
   REQUIRE(i3(-1.0) == 0.0);
   REQUIRE(i3(0.0) == 0.0);
   REQUIRE(i3(0.25) == 0.25);
   REQUIRE(i3(0.50) == 0.25 * 1.5);
   REQUIRE(i3(0.75) == Approx(1.25 / 2.0));
   REQUIRE(i3(1.00) == 0.5);
   REQUIRE(i3(2.00) == 0.0);
   REQUIRE(i4(-1.0) == 0.0);
   REQUIRE(i4(0.0) == 0.0);
   REQUIRE(i4(0.25) == 0.25);
   REQUIRE(i4(0.50) == 0.25 * 1.5);
   REQUIRE(i4(0.75) == Approx(1.25 / 2.0));
   REQUIRE(i4(1.00) == 0.5);
   REQUIRE(i4(2.00) == 0.0);
}

TEST_CASE("Verify quotient of interpolants.", "[interpolant]")
{
   ilist<double, double> list1 = {{0.00, 0.00}, {0.50, 0.25}, {1.00, 1.00}};
   ilist<double, double> list2 = {{0.25, 2.00}, {0.75, 1.00}, {1.25, 0.50}};
   interpolantd i1(list1);
   interpolantd i2(list2);
   interpolantd i3(i1 / i2);
   REQUIRE(i3(-1.0) == 0.0);
   REQUIRE(i3(0.0) == 0.0);
   REQUIRE(i3(0.25) == (0.25 / 2.0) / 2.0);
   REQUIRE(i3(0.50) == 0.25 / 1.5);
   REQUIRE(i3(0.75) == Approx(1.25 / 2.0));
   REQUIRE(i3(1.00) == 1.0 / 0.75);
   REQUIRE(i3(2.00) == 1.0 / 0.5);
}

TEST_CASE("Verify file input.", "[interpolant]")
{
   REQUIRE_THROWS(interpolantd j("nonexistent-file"));
   REQUIRE_THROWS(interpolantd k("interpolant_test-badinput.txt"));
   REQUIRE_THROWS(interpolantd m("interpolant_test-badinput2.txt"));
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
   interpolantd n("interpolant_test-shortinput.txt");
   cout << "n.points().size()=" << n.points().size() << endl;
   REQUIRE(n.integral() == 0.0);
}

TEST_CASE("Verify file input and interoperation with units.", "[interpolant]")
{
   interpolant<length, num::time> i("interpolant_test.txt", nm, s);
   interpolant<dyndim, dyndim> j("interpolant_test.txt", nm, s);

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

   REQUIRE(j(-1 * nm) == 0.0 * s);
   REQUIRE(j(299 * nm) == 0.0 * s);
   REQUIRE(j(300 * nm) == 0.0 * s);
   REQUIRE((j(350 * nm) / s).number() == Approx(40.0));
   REQUIRE(j(400 * nm) == 80.0 * s);
   REQUIRE(j(450 * nm) == 80.0 * s);
   REQUIRE(j(500 * nm) == 80.0 * s);
   REQUIRE(j(900 * nm) == 80.0 * s);
   REQUIRE((j(950 * nm) / s).number() == Approx(40.0));
   REQUIRE(j(1000 * nm) == 0.0 * s);
   REQUIRE(j(1100 * nm) == 0.0 * s);
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

   interpolant<length, double> i(list1);
   interpolant<dyndim, dyndim> j(list2);

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

   REQUIRE(j(-1 * nm) == 0.0 * s);
   REQUIRE(j(299 * nm) == 0.0 * s);
   REQUIRE(j(300 * nm) == 0.0 * s);
   REQUIRE((j(350 * nm) / s).number() == Approx(7.5));
   REQUIRE(j(400 * nm) == 15.0 * s);
   REQUIRE((j(450 * nm) / s).number() == Approx(0.5 * (15.0 + 80.0)));
   REQUIRE(j(500 * nm) == 80.0 * s);
   REQUIRE(j(600 * nm) == 80.0 * s);
   REQUIRE(j(900 * nm) == 80.0 * s);
   REQUIRE((j(950 * nm) / s).number() == Approx(40.0));
   REQUIRE(j(1000 * nm) == 0.0 * s);
   REQUIRE(j(1100 * nm) == 0.0 * s);
}

TEST_CASE("Verify scalar multiplication of interpolant on right.",
          "[interpolant]")
{
   interpolant<length, num::time> i1("interpolant_test.txt", nm, s);
   interpolant<dyndim, dyndim> i2("interpolant_test.txt", nm, s);
   interpolant<length, length> j1 = i1 * (2 * m / s);
   interpolant<dyndim, dyndim> j2 = i2 * (2 * m / s);

   REQUIRE(j1(-1 * nm) == 0.0 * m);
   REQUIRE(j1(299 * nm) == 0.0 * m);
   REQUIRE(j1(300 * nm) == 0.0 * m);
   REQUIRE(j1(350 * nm) / m == Approx(80.0));
   REQUIRE(j1(400 * nm) == 160.0 * m);
   REQUIRE(j1(450 * nm) == 160.0 * m);
   REQUIRE(j1(500 * nm) == 160.0 * m);
   REQUIRE(j1(900 * nm) == 160.0 * m);
   REQUIRE(j1(950 * nm) / m == Approx(80.0));
   REQUIRE(j1(1000 * nm) == 0.0 * m);
   REQUIRE(j1(1100 * nm) == 0.0 * m);

   REQUIRE(j2(-1 * nm) == 0.0 * m);
   REQUIRE(j2(299 * nm) == 0.0 * m);
   REQUIRE(j2(300 * nm) == 0.0 * m);
   REQUIRE((j2(350 * nm) / m).number() == Approx(80.0));
   REQUIRE(j2(400 * nm) == 160.0 * m);
   REQUIRE(j2(450 * nm) == 160.0 * m);
   REQUIRE(j2(500 * nm) == 160.0 * m);
   REQUIRE(j2(900 * nm) == 160.0 * m);
   REQUIRE((j2(950 * nm) / m).number() == Approx(80.0));
   REQUIRE(j2(1000 * nm) == 0.0 * m);
   REQUIRE(j2(1100 * nm) == 0.0 * m);

   interpolant<length, num::time> k1 = i1 * 2.0;
   interpolant<dyndim, dyndim> k2 = i2 * 2.0;

   REQUIRE(k1(-1 * nm) == 0.0 * s);
   REQUIRE(k1(299 * nm) == 0.0 * s);
   REQUIRE(k1(300 * nm) == 0.0 * s);
   REQUIRE(k1(350 * nm) / s == Approx(80.0));
   REQUIRE(k1(400 * nm) == 160.0 * s);
   REQUIRE(k1(450 * nm) == 160.0 * s);
   REQUIRE(k1(500 * nm) == 160.0 * s);
   REQUIRE(k1(900 * nm) == 160.0 * s);
   REQUIRE(k1(950 * nm) / s == Approx(80.0));
   REQUIRE(k1(1000 * nm) == 0.0 * s);
   REQUIRE(k1(1100 * nm) == 0.0 * s);

   REQUIRE(k2(-1 * nm) == 0.0 * s);
   REQUIRE(k2(299 * nm) == 0.0 * s);
   REQUIRE(k2(300 * nm) == 0.0 * s);
   REQUIRE((k2(350 * nm) / s).number() == Approx(80.0));
   REQUIRE(k2(400 * nm) == 160.0 * s);
   REQUIRE(k2(450 * nm) == 160.0 * s);
   REQUIRE(k2(500 * nm) == 160.0 * s);
   REQUIRE(k2(900 * nm) == 160.0 * s);
   REQUIRE((k2(950 * nm) / s).number() == Approx(80.0));
   REQUIRE(k2(1000 * nm) == 0.0 * s);
   REQUIRE(k2(1100 * nm) == 0.0 * s);

   i1 *= 2.0;
   i2 *= 2.0;

   REQUIRE(i1(-1 * nm) == 0.0 * s);
   REQUIRE(i1(299 * nm) == 0.0 * s);
   REQUIRE(i1(300 * nm) == 0.0 * s);
   REQUIRE(i1(350 * nm) / s == Approx(80.0));
   REQUIRE(i1(400 * nm) == 160.0 * s);
   REQUIRE(i1(450 * nm) == 160.0 * s);
   REQUIRE(i1(500 * nm) == 160.0 * s);
   REQUIRE(i1(900 * nm) == 160.0 * s);
   REQUIRE(i1(950 * nm) / s == Approx(80.0));
   REQUIRE(i1(1000 * nm) == 0.0 * s);
   REQUIRE(i1(1100 * nm) == 0.0 * s);

   REQUIRE(i2(-1 * nm) == 0.0 * s);
   REQUIRE(i2(299 * nm) == 0.0 * s);
   REQUIRE(i2(300 * nm) == 0.0 * s);
   REQUIRE((i2(350 * nm) / s).number() == Approx(80.0));
   REQUIRE(i2(400 * nm) == 160.0 * s);
   REQUIRE(i2(450 * nm) == 160.0 * s);
   REQUIRE(i2(500 * nm) == 160.0 * s);
   REQUIRE(i2(900 * nm) == 160.0 * s);
   REQUIRE((i2(950 * nm) / s).number() == Approx(80.0));
   REQUIRE(i2(1000 * nm) == 0.0 * s);
   REQUIRE(i2(1100 * nm) == 0.0 * s);
}

TEST_CASE("Verify scalar multiplication of interpolant on left.",
          "[interpolant]")
{
   interpolant<length, num::time> i1("interpolant_test.txt", nm, s);
   interpolant<dyndim, dyndim> i2("interpolant_test.txt", nm, s);
   interpolant<length, length> j1 = (2 * m / s) * i1;
   interpolant<dyndim, dyndim> j2 = (2 * m / s) * i2;

   REQUIRE(j1(-1 * nm) == 0.0 * m);
   REQUIRE(j1(299 * nm) == 0.0 * m);
   REQUIRE(j1(300 * nm) == 0.0 * m);
   REQUIRE(j1(350 * nm) / m == Approx(80.0));
   REQUIRE(j1(400 * nm) == 160.0 * m);
   REQUIRE(j1(450 * nm) == 160.0 * m);
   REQUIRE(j1(500 * nm) == 160.0 * m);
   REQUIRE(j1(900 * nm) == 160.0 * m);
   REQUIRE(j1(950 * nm) / m == Approx(80.0));
   REQUIRE(j1(1000 * nm) == 0.0 * m);
   REQUIRE(j1(1100 * nm) == 0.0 * m);

   REQUIRE(j2(-1 * nm) == 0.0 * m);
   REQUIRE(j2(299 * nm) == 0.0 * m);
   REQUIRE(j2(300 * nm) == 0.0 * m);
   REQUIRE((j2(350 * nm) / m).number() == Approx(80.0));
   REQUIRE(j2(400 * nm) == 160.0 * m);
   REQUIRE(j2(450 * nm) == 160.0 * m);
   REQUIRE(j2(500 * nm) == 160.0 * m);
   REQUIRE(j2(900 * nm) == 160.0 * m);
   REQUIRE((j2(950 * nm) / m).number() == Approx(80.0));
   REQUIRE(j2(1000 * nm) == 0.0 * m);
   REQUIRE(j2(1100 * nm) == 0.0 * m);

   interpolant<length, num::time> k1 = 2.0 * i1;
   interpolant<dyndim, dyndim> k2 = 2.0 * i2;

   REQUIRE(k1(-1 * nm) == 0.0 * s);
   REQUIRE(k1(299 * nm) == 0.0 * s);
   REQUIRE(k1(300 * nm) == 0.0 * s);
   REQUIRE(k1(350 * nm) / s == Approx(80.0));
   REQUIRE(k1(400 * nm) == 160.0 * s);
   REQUIRE(k1(450 * nm) == 160.0 * s);
   REQUIRE(k1(500 * nm) == 160.0 * s);
   REQUIRE(k1(900 * nm) == 160.0 * s);
   REQUIRE(k1(950 * nm) / s == Approx(80.0));
   REQUIRE(k1(1000 * nm) == 0.0 * s);
   REQUIRE(k1(1100 * nm) == 0.0 * s);

   REQUIRE(k2(-1 * nm) == 0.0 * s);
   REQUIRE(k2(299 * nm) == 0.0 * s);
   REQUIRE(k2(300 * nm) == 0.0 * s);
   REQUIRE((k2(350 * nm) / s).number() == Approx(80.0));
   REQUIRE(k2(400 * nm) == 160.0 * s);
   REQUIRE(k2(450 * nm) == 160.0 * s);
   REQUIRE(k2(500 * nm) == 160.0 * s);
   REQUIRE(k2(900 * nm) == 160.0 * s);
   REQUIRE((k2(950 * nm) / s).number() == Approx(80.0));
   REQUIRE(k2(1000 * nm) == 0.0 * s);
   REQUIRE(k2(1100 * nm) == 0.0 * s);
}

TEST_CASE("Verify scalar division of interpolant.", "[interpolant]")
{
   interpolant<length, num::time> i1("interpolant_test.txt", nm, s);
   interpolant<dyndim, dyndim> i2("interpolant_test.txt", nm, s);
   interpolant<length, double> j1 = i1 / (2 * s);
   interpolant<dyndim, dyndim> j2 = i2 / (2 * s);

   REQUIRE(j1(-1 * nm) == 0.0);
   REQUIRE(j1(299 * nm) == 0.0);
   REQUIRE(j1(300 * nm) == 0.0);
   REQUIRE(j1(350 * nm) == Approx(20.0));
   REQUIRE(j1(400 * nm) == 40.0);
   REQUIRE(j1(450 * nm) == 40.0);
   REQUIRE(j1(500 * nm) == 40.0);
   REQUIRE(j1(900 * nm) == 40.0);
   REQUIRE(j1(950 * nm) == Approx(20.0));
   REQUIRE(j1(1000 * nm) == 0.0);
   REQUIRE(j1(1100 * nm) == 0.0);

   REQUIRE((j2(-1 * nm)).number() == 0.0);
   REQUIRE((j2(299 * nm)).number() == 0.0);
   REQUIRE((j2(300 * nm)).number() == 0.0);
   REQUIRE((j2(350 * nm)).number() == Approx(20.0));
   REQUIRE((j2(400 * nm)).number() == 40.0);
   REQUIRE((j2(450 * nm)).number() == 40.0);
   REQUIRE((j2(500 * nm)).number() == 40.0);
   REQUIRE((j2(900 * nm)).number() == 40.0);
   REQUIRE((j2(950 * nm)).number() == Approx(20.0));
   REQUIRE((j2(1000 * nm)).number() == 0.0);
   REQUIRE((j2(1100 * nm)).number() == 0.0);

   i1 /= 2.0;
   i2 /= 2.0;

   REQUIRE(i1(-1 * nm) == 0.0 * s);
   REQUIRE(i1(299 * nm) == 0.0 * s);
   REQUIRE(i1(300 * nm) == 0.0 * s);
   REQUIRE(i1(350 * nm) / s == Approx(20.0));
   REQUIRE(i1(400 * nm) == 40.0 * s);
   REQUIRE(i1(450 * nm) == 40.0 * s);
   REQUIRE(i1(500 * nm) == 40.0 * s);
   REQUIRE(i1(900 * nm) == 40.0 * s);
   REQUIRE(i1(950 * nm) / s == Approx(20.0));
   REQUIRE(i1(1000 * nm) == 0.0 * s);
   REQUIRE(i1(1100 * nm) == 0.0 * s);

   REQUIRE(i2(-1 * nm) == 0.0 * s);
   REQUIRE(i2(299 * nm) == 0.0 * s);
   REQUIRE(i2(300 * nm) == 0.0 * s);
   REQUIRE((i2(350 * nm) / s).number() == Approx(20.0));
   REQUIRE(i2(400 * nm) == 40.0 * s);
   REQUIRE(i2(450 * nm) == 40.0 * s);
   REQUIRE(i2(500 * nm) == 40.0 * s);
   REQUIRE(i2(900 * nm) == 40.0 * s);
   REQUIRE((i2(950 * nm) / s).number() == Approx(20.0));
   REQUIRE(i2(1000 * nm) == 0.0 * s);
   REQUIRE(i2(1100 * nm) == 0.0 * s);
}

TEST_CASE(
      "Verify that multiplication of interpolant by empty interpolant is "
      "empty",
      "[interpolant]")
{
   interpolant<length, num::time> i1("interpolant_test.txt", nm, s);
   interpolant<dyndim, dyndim> i2("interpolant_test.txt", nm, s);
   interpolant<length, double> j1;
   interpolant<dyndim, dyndim> j2;
   REQUIRE_THROWS((i1 * j1)(0 * cm));
   REQUIRE_THROWS((i2 * j2)(0 * cm));
}

double my_erf(double x) { return erf(x); }

TEST_CASE("Verify that interpolant of function produces right integral.",
          "[interpolant]")
{
   function<area(length)> square = [](length x) { return x * x; };
   volume const i = integral(square, 0 * cm, 1 * cm);

   interpolant<length, area> const j1(square, 0 * cm, 1 * cm, 1.0E-17);
   interpolant<dyndim, dyndim> const j2(square, 0 * cm, 1 * cm, 1.0E-17);

   REQUIRE(j1.integral() / i == Approx(1.0));
   REQUIRE(j1.integral(0 * cm, 1 * cm) / i == Approx(1.0));

   REQUIRE((j2.integral() / i).number() == Approx(1.0));
   REQUIRE((j2.integral(0 * cm, 1 * cm) / i).number() == Approx(1.0));

   interpolant<length, area> const k1(square, 1 * cm, 0 * cm);
   interpolant<dyndim, dyndim> const k2(square, 1 * cm, 0 * cm);

   REQUIRE(k1.integral(0 * cm, 1 * cm) / i == Approx(1.0));
   REQUIRE((k2.integral(0 * cm, 1 * cm) / i).number() == Approx(1.0));

   REQUIRE_THROWS(
         (interpolant<length, area>(square, 0 * cm, 1 * cm, -1.0E-06)));

   REQUIRE_THROWS(
         (interpolant<dyndim, dyndim>(square, 0 * cm, 1 * cm, -1.0E-06)));

   double tol = 1.0E-04;
   interpolantd const e(erf, -1.0, 2.0, tol);
   REQUIRE(e.integral(-1.0, 2.0) / integral(my_erf, -1.0, 2.0, tol) ==
           Approx(1.0).epsilon(tol));

   tol = 1.0E-03;
   function<double(double)> g = [](double x) { return exp(-0.5 * x * x); };
   interpolantd const ig(g, -5.0, +5.0, tol);
   REQUIRE(ig.integral() == Approx(sqrt(2.0 * M_PI)).epsilon(tol));
   ofstream ofg("gaussian.dat");
   ofg << ig.points();
}

