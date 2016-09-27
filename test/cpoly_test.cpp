
// Copyright 2016  Thomas E. Vaughan
//
// This software is distributable under the terms of the GNU LGPL, Version 3 or
// later.

#include <sstream> // for ostringstream

#include "catch.hpp"
#include "cpoly.hpp"
#include "units.hpp"

using namespace num;

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

