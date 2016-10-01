
// Copyright 2016  Thomas E. Vaughan
//
// This software is distributable under the terms of the GNU LGPL, Version 3 or
// later.

#include <sstream> // for ostringstream

#include "catch.hpp"
#include "cpoly.hpp"
#include "units.hpp"

using namespace std;
using namespace num;
using namespace num::u;

TEST_CASE("Verify default construction.", "[cpoly]")
{
   cpoly<4> cp1;
   REQUIRE(cp1.num_coefs() == 5);
   REQUIRE(cp1.coef<0>() == 0);
   REQUIRE(cp1.coef<1>() == 0);
   REQUIRE(cp1.coef<2>() == 0);
   REQUIRE(cp1.coef<3>() == 0);
   REQUIRE(cp1.coef<4>() == 0);

   cpoly<2, dyndim, dyndim> cp2;
   REQUIRE(cp2.num_coefs() == 3);
   REQUIRE(cp2.coef<0>().number() == 0);
   REQUIRE(cp2.coef<1>().number() == 0);
   REQUIRE(cp2.coef<2>().number() == 0);

   cpoly<2, num::time, length> cp3;
   REQUIRE(cp3.num_coefs() == 3);
   REQUIRE(cp3.coef<0>() == 0 * m);
   REQUIRE(cp3.coef<1>() == 0 * m / s);
   REQUIRE(cp3.coef<2>() == 0 * m / s / s);
}

TEST_CASE("Verify construction from array.", "[cpoly]")
{
   array<dyndim, 3> a1;
   vector<dyndim> v1(3);
   v1[0] = a1[0] = 1 * m;
   v1[1] = a1[1] = 1 * m / s;
   v1[2] = a1[2] = 0.5 * m / s / s;
   cpoly<2, dyndim, dyndim> cp1(a1);
   cpoly<2, dyndim, dyndim> cp2(v1);
   REQUIRE(cp1.coef<0>() == a1[0]);
   REQUIRE(cp1.coef<1>() == a1[1]);
   REQUIRE(cp1.coef<2>() == a1[2]);
   REQUIRE(cp1.coef<0>() == v1[0]);
   REQUIRE(cp1.coef<1>() == v1[1]);
   REQUIRE(cp1.coef<2>() == v1[2]);
}

TEST_CASE("Verify evaluation of polynomial.", "[cpoly]")
{
   array<dyndim, 3> a1;
   a1[0] = 1 * m;
   a1[1] = 1 * m / s;
   a1[2] = 0.5 * m / s / s;
   cpoly<2, dyndim, dyndim> cp1(a1);
   REQUIRE(cp1(0 * s) == 1.0 * m);
   REQUIRE(cp1(1 * s) == 2.5 * m);
   REQUIRE(cp1(2 * s) == 5.0 * m);
}

TEST_CASE("Verify derivative of polynomial.", "[cpoly]")
{
   array<dyndim, 3> a1;
   a1[0] = 1 * m;
   a1[1] = 1 * m / s;
   a1[2] = 0.5 * m / s / s;
   cpoly<2, dyndim, dyndim> cp1(a1);
   auto cp2 = cp1.derivative();
   REQUIRE(cp2.num_coefs() == 2);
   REQUIRE(cp2.coef<0>() == 1.0 * m / s);
   REQUIRE(cp2.coef<1>() == 1.0 * m / s / s);
}

TEST_CASE("Verify integral of polynomial.", "[cpoly]")
{
   array<dyndim, 3> a1;
   a1[0] = 1 * m;
   a1[1] = 1 * m / s;
   a1[2] = 0.5 * m / s / s;
   cpoly<2, dyndim, dyndim> cp1(a1);
   auto cp2 = cp1.derivative();
   auto cp3 = cp2.integral(1 * s); // Start integrating at 1 sec.
   REQUIRE(cp3.num_coefs() == 3);
   REQUIRE(cp3.coef<0>() == -1.5 * m);
   REQUIRE(cp3.coef<1>() == +1.0 * m / s);
   REQUIRE(cp3.coef<2>() == +0.5 * m / s / s);
}

TEST_CASE("Verify derivative down to constant and integral back to linear.",
          "[cpoly]")
{
   array<dyndim, 3> a1;
   a1[0] = 1 * m;
   a1[1] = 1 * m / s;
   a1[2] = 0.5 * m / s / s;
   cpoly<2, dyndim, dyndim> cp1(a1); // quadratic
   auto cp2 = cp1.derivative();      // linear
   auto cp3 = cp2.derivative();      // constant
   REQUIRE(acceleration(cp3) == 1 * m / s / s);
   REQUIRE(cp3(1 * s) == 1 * m / s / s);
   REQUIRE(cp3(2 * s) == 1 * m / s / s);
   // Change constant.
   cp3 = 2 * m / s / s;
   auto cp4 = cp3.integral(0.5 * s);
   REQUIRE(cp4.num_coefs() == 2);
   REQUIRE(cp4.coef<0>() == -1 * m / s);
   REQUIRE(cp4.coef<1>() == +2 * m / s / s);
}

