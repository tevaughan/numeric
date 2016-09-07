
// Copyright 2016  Thomas E. Vaughan
//
// This software is distributable under the terms of the GNU LGPL, Version 3 or
// later.

#include <sstream> // for ostringstream

#include "catch.hpp"
#include "units.hpp"

using namespace num;
using namespace std;

TEST_CASE("Verify basic operation of units.", "[units]")
{
   length const len1 = meters(3.0);
   length const len2 = 6.0 * m;
   REQUIRE(len1 == 3.0 * m);
   REQUIRE(len2 == meters(6.0));
   REQUIRE(len1 / len2 == Approx(0.5));
}

TEST_CASE("Verify string conversion.", "[units]")
{
   newtons const f1(2.0);
   force const f2 = 3.0 * N;
   ostringstream oss1, oss2;
   oss1 << f1;
   oss2 << f2;
   REQUIRE(oss1.str() == "[2 N]");
   REQUIRE(oss2.str() == "[3 kg m s^-2]");

   microkelvins const t1(2.5);
   temperature const t2 = 2.5 * muK;
   ostringstream oss3, oss4;
   oss3 << t1;
   oss4 << t2;
   REQUIRE(oss3.str() == "[2.5 muK]");
   REQUIRE(oss4.str() == "[2.5e-06 K]");

   dimval<1,1,1,1,1> const foo = s * m * kg * C * K;
   ostringstream oss5, oss6;
   oss5 << foo;
   oss6 << foo.pow<-1>();
   REQUIRE(oss5.str() == "[1 kg m s C K]");
   REQUIRE(oss6.str() == "[1 kg^-1 m^-1 s^-1 C^-1 K^-1]");
}

TEST_CASE("Verify default construction of zero value.", "[units]")
{
   REQUIRE(length() == 0 * m);
}

TEST_CASE("Verify addition.", "[units]") { REQUIRE(m + km == 1001 * m); }

TEST_CASE("Verify subtraction.", "[units]")
{
   REQUIRE(10 * K - 1 * K == 9 * K);
}

TEST_CASE("Verify additive assignment.", "[units]")
{
   speed v1 = 2 * m / s;
   speed const v2 = 3 * m / s;
   REQUIRE((v1 += v2) == 5 * m / s);
}

TEST_CASE("Verify subtractive assignment.", "[units]")
{
   charge q1 = 3 * C;
   charge const q2 = 1 * C;
   REQUIRE((q1 -= q2) == 2 * C);
}

TEST_CASE("Verify multiplication of dimvals.", "[units]")
{
   force const f = 10 * N;
   length const x = 10 * m;
   REQUIRE(f * x == 100 * J);
}

TEST_CASE("Verify multiplication of dimval by number on right.", "[units]")
{
   num::time const t = 2 * s;
   REQUIRE(t * 2 == 4 * s);
}

TEST_CASE("Verify multiplication of dimval by number on left.", "[units]")
{
   power const p = 2 * W;
   REQUIRE(2 * p == 4 * W);
}

TEST_CASE("Verify division resulting in number.", "[units]")
{
   length const x1 = 1 * km;
   length const x2 = 1 * m;
   REQUIRE(x1 / x2 == 1000);
}

TEST_CASE("Verify division of dimvals.", "[units]")
{
   energy const e = 6 * J;
   num::time const t = 3 * s;
   REQUIRE(e / t == 2 * W);
}

TEST_CASE("Verify division of number by dimval.", "[units]")
{
   double const n = 3;
   frequency const f = 1 * Hz;
   REQUIRE(n / f == 3 * s);
}

TEST_CASE("Verify division of dimval by number.", "[units]")
{
   double const n = 3;
   frequency const f = 6 * Hz;
   REQUIRE(f / n == 2 * Hz);
}

TEST_CASE("Verify multiplicative assignment.", "[units]")
{
   double const n = 3;
   area a = m * m;
   REQUIRE((a *= n) == 3 * m * m);
}

TEST_CASE("Verify divisive assignment.", "[units]")
{
   double const n = 3;
   energy e = 3 * J;
   REQUIRE((e /= n) == 1 * J);
}

TEST_CASE("Verify comparison.", "[units]")
{
   mass const m1 = 1 * kg;
   mass const m2 = 2 * kg;
   REQUIRE(m1 < m2);
   REQUIRE(m2 > m1);
   REQUIRE(m1 <= m1);
   REQUIRE(m1 >= m1);
   REQUIRE(m1 == m1);
   REQUIRE(m1 != m2);
}

TEST_CASE("Verify integer power.", "[units]")
{
   length const x = 2 * km;
   REQUIRE(x.pow<2>() == 4 * km * km);
   REQUIRE(x.pow<-2>() * (km * km) == Approx(0.25));
}

TEST_CASE("Verify integer root.", "[units]")
{
   volume const v = 8 * km * km * km;
   REQUIRE(v.root<3>() / km == Approx(2.0));
   REQUIRE(sqrt(v * 2 * km) / (km * km) == Approx(4.0));
}

TEST_CASE("Verify absolute value.", "[units]")
{
   length const x = -1 * m;
   REQUIRE(fabs(x) == 1 * m);
}

