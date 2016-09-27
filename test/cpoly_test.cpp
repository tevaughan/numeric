
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
   REQUIRE(cp1.c().size() == 5);
   REQUIRE(cp1.c()[0] == 0);
   REQUIRE(cp1.c()[1] == 0);
   REQUIRE(cp1.c()[2] == 0);
   REQUIRE(cp1.c()[3] == 0);
   REQUIRE(cp1.c()[4] == 0);

   cpoly<2, dyndim, dyndim> cp2;
   REQUIRE(cp2.c().size() == 3);
   REQUIRE(cp2.c()[0].number() == 0);
   REQUIRE(cp2.c()[1].number() == 0);
   REQUIRE(cp2.c()[2].number() == 0);

   cpoly<2, dyndim, length> cp3;
   REQUIRE(cp3.c().size() == 3);
   REQUIRE(cp3.c()[0].number() == 0);
   REQUIRE(cp3.c()[1].number() == 0);
   REQUIRE(cp3.c()[2].number() == 0);
}

TEST_CASE("Verify construction from array.", "[cpoly]")
{
   array<dyndim, 3> a1;
   vector<dyndim> v1(3);
   v1[0] = a1[0] = 1 * m;
   v1[1] = a1[1] = 1 * m / s;
   v1[2] = a1[2] = 0.5 * m / s / s;
   cpoly<2, dyndim, length> cp1(a1);
   cpoly<2, dyndim, length> cp2(v1);
   REQUIRE(cp1.c()[0] == a1[0]);
   REQUIRE(cp1.c()[1] == a1[1]);
   REQUIRE(cp1.c()[2] == a1[2]);
   REQUIRE(cp1.c()[0] == v1[0]);
   REQUIRE(cp1.c()[1] == v1[1]);
   REQUIRE(cp1.c()[2] == v1[2]);
}

TEST_CASE("Verify evaltuation of polynomial.", "[cpoly]")
{
   array<dyndim, 3> a1;
   a1[0] = 1 * m;
   a1[1] = 1 * m / s;
   a1[2] = 0.5 * m / s / s;
   cpoly<2, dyndim, length> cp1(a1);
   REQUIRE(cp1(0 * s) == 1.0 * m);
   REQUIRE(cp1(1 * s) == 2.5 * m);
   REQUIRE(cp1(2 * s) == 5.0 * m);
}

