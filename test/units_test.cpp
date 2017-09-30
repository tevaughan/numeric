
// Copyright 2016-2017  Thomas E. Vaughan
//
// This software is distributable under the terms of the GNU LGPL, Version 3 or
// later.

#include <sstream> // for ostringstream

#include "catch.hpp"
#include "units.hpp"

using namespace num;
using namespace num::u;
using namespace std;

TEST_CASE("Verify basic operation of units.", "[units]")
{
   length const len1 = meters(3.0);
   length const len2 = 6.0 * m;
   dyndim const len3 = len1;
   dyndim const len4 = 6.0 * m;
   REQUIRE(len1 == 3.0 * m);
   REQUIRE(len2 == meters(6.0));
   REQUIRE(len1 / len2 == Approx(0.5));
   REQUIRE(len3 == len1);
   REQUIRE(len4 == meters(6.0));
}

TEST_CASE("Verify string conversion.", "[units]")
{
   newtons const f1(2.0);
   force const f2 = 3.0 * N;
   dyndim const f3 = newtons(2.0);
   dyndim const f4 = 3.0 * N;
   ostringstream oss1, oss2, ossf3, ossf4;
   oss1 << f1;
   oss2 << f2;
   ossf3 << f3;
   ossf4 << f4;
   REQUIRE(oss1.str() == "[2 N]");
   REQUIRE(oss2.str() == "[3 kg m s^-2]");
   REQUIRE(ossf3.str() == "[2 kg m s^-2]");
   REQUIRE(ossf4.str() == "[3 kg m s^-2]");

   microkelvins const t1(2.5);
   temperature const t2 = 2.5 * muK;
   dyndim const t3 = microkelvins(2.5);
   dyndim const t4 = 2.5 * muK;
   ostringstream oss3, oss4, osst3, osst4;
   oss3 << t1;
   oss4 << t2;
   osst3 << t3;
   osst4 << t4;
   REQUIRE(oss3.str() == "[2.5 muK]");
   REQUIRE(oss4.str() == "[2.5e-06 K]");
   REQUIRE(osst3.str() == "[2.5e-06 K]");
   REQUIRE(osst4.str() == "[2.5e-06 K]");

   statdim<1, 1, 1, 1, 1> const foo = s * m * kg * C * K;
   dyndim const bar = s * m * kg * C * K;
   ostringstream oss5, oss6, oss7, oss8;
   oss5 << foo;
   oss6 << foo.pow<-1>();
   oss7 << bar;
   oss8 << bar.pow(-1);
   REQUIRE(oss5.str() == "[1 kg m s C K]");
   REQUIRE(oss6.str() == "[1 kg^-1 m^-1 s^-1 C^-1 K^-1]");
   REQUIRE(oss7.str() == "[1 kg m s C K]");
   REQUIRE(oss8.str() == "[1 kg^-1 m^-1 s^-1 C^-1 K^-1]");
}

TEST_CASE("Verify default construction of zero value.", "[units]")
{
   REQUIRE(length() == 0 * m);
}

TEST_CASE("Verify addition.", "[units]") { REQUIRE(m + km == 1001 * m); }

TEST_CASE("Verify subtraction.", "[units]")
{
   REQUIRE(10 * K - 1 * K == 9 * K); // statdim - statdim
   dyndim const x = 10 * K;
   dyndim const y = 1 * K;
   REQUIRE(x - y == 9 * K);      // dyndim - dyndim
   REQUIRE(10 * K - x == 0 * K); // statdim - dyndim
   REQUIRE(y - 1 * K == 0 * K);  // dyndim - statdim
}

TEST_CASE("Verify additive assignment.", "[units]")
{
   speed v1 = 2 * m / s;
   speed const v2 = 3 * m / s;
   REQUIRE((v1 += v2) == 5 * m / s); // statdim += statdim
   dyndim v3 = v1;
   dyndim const v4 = v2;
   REQUIRE((v1 += v4) == 8 * m / s);  // statdim += dyndim
   REQUIRE((v3 += v2) == 8 * m / s);  // dyndim += statdim
   REQUIRE((v3 += v4) == 11 * m / s); // dyndim += dyndim
}

TEST_CASE("Verify subtractive assignment.", "[units]")
{
   charge q1 = 3 * C;
   charge const q2 = 1 * C;
   REQUIRE((q1 -= q2) == 2 * C); // statdim -= statdim
   dyndim q3 = q1;
   dyndim const q4 = q2;
   REQUIRE((q1 -= q4) == 1 * C); // statdim -= dyndim
   REQUIRE((q3 -= q2) == 1 * C); // dyndim -= statdim
   REQUIRE((q3 -= q4) == 0 * C); // dyndim -= dyndim
}

TEST_CASE("Verify multiplication of dimvals.", "[units]")
{
   force const f1 = 10 * N;
   dyndim const f2 = 10 * N;
   length const x1 = 10 * m;
   dyndim const x2 = 10 * m;
   REQUIRE(f1 * x1 == 100 * J); // statdim * statdim
   REQUIRE(f1 * x2 == 100 * J); // statdim * dyndim
   REQUIRE(f2 * x1 == 100 * J); // dyndim * statdim
   REQUIRE(f2 * x2 == 100 * J); // dyndim * dyndim
}

TEST_CASE("Verify multiplication of dimval by number on right.", "[units]")
{
   num::time const t1 = 2 * s;
   dyndim const t2 = 2 * s;
   REQUIRE(t1 * 2 == 4 * s); // statdim * number
   REQUIRE(t2 * 2 == 4 * s); // dyndim * number
}

TEST_CASE("Verify multiplication of dimval by number on left.", "[units]")
{
   power const p1 = 2 * W;
   dyndim const p2 = 2 * W;
   REQUIRE(2 * p1 == 4 * W); // number * statdim
   REQUIRE(2 * p2 == 4 * W); // number * dyndim
}

TEST_CASE("Verify division resulting in number.", "[units]")
{
   length const x1 = 1 * km;
   dyndim const x2 = 1 * km;
   length const y1 = 1 * m;
   dyndim const y2 = 1 * m;
   REQUIRE(x1 / y1 == 1000);            // statdim-length / statdim-length
   REQUIRE((x1 / y2).number() == 1000); // statdim-length / dyndim-length
   REQUIRE((x2 / y1).number() == 1000); // dyndim-length / statdim-length
   REQUIRE((x2 / y2).number() == 1000); // dyndim-length / dyndim-length
}

TEST_CASE("Verify division of dimvals.", "[units]")
{
   energy const e1 = 6 * J;
   dyndim const e2 = 6 * J;
   num::time const t1 = 3 * s;
   dyndim const t2 = 3 * s;
   REQUIRE(e1 / t1 == 2 * W); // statdim / statdim
   REQUIRE(e1 / t2 == 2 * W); // statdim / dyndim
   REQUIRE(e2 / t1 == 2 * W); // dyndim / statdim
   REQUIRE(e2 / t2 == 2 * W); // dyndim / dyndim
}

TEST_CASE("Verify division of number by dimval.", "[units]")
{
   double const n = 3;
   frequency const f1 = 1 * Hz;
   dyndim const f2 = 1 * Hz;
   REQUIRE(n / f1 == 3 * s); // number / statdim
   REQUIRE(n / f2 == 3 * s); // number / dyndim
}

TEST_CASE("Verify division of dimval by number.", "[units]")
{
   double const n = 3;
   frequency const f1 = 6 * Hz;
   dyndim const f2 = 6 * Hz;
   REQUIRE(f1 / n == 2 * Hz); // statdim / number
   REQUIRE(f2 / n == 2 * Hz); // dyndim / number
}

TEST_CASE("Verify multiplicative assignment.", "[units]")
{
   double const n = 3;
   area a1 = m * m;
   dyndim a2 = m * m;
   REQUIRE((a1 *= n) == 3 * m * m);
   REQUIRE((a2 *= n) == 3 * m * m);
}

TEST_CASE("Verify divisive assignment.", "[units]")
{
   double const n = 3;
   energy e1 = 3 * J;
   dyndim e2 = 3 * J;
   REQUIRE((e1 /= n) == 1 * J);
   REQUIRE((e2 /= n) == 1 * J);
}

TEST_CASE("Verify comparison.", "[units]")
{
   mass const m1 = 1 * kg;
   dyndim const m2 = 1 * kg;
   mass const n1 = 2 * kg;
   dyndim const n2 = 2 * kg;

   REQUIRE(m1 < n1);
   REQUIRE(m2 < n1);
   REQUIRE(m1 < n2);
   REQUIRE(m2 < n2);

   REQUIRE(n1 > m1);
   REQUIRE(n2 > m1);
   REQUIRE(n1 > m2);
   REQUIRE(n2 > m2);

   REQUIRE(m1 <= m1);
   REQUIRE(m2 <= m1);
   REQUIRE(m1 <= m2);
   REQUIRE(m2 <= m2);

   REQUIRE(m1 >= m1);
   REQUIRE(m2 >= m1);
   REQUIRE(m1 >= m2);
   REQUIRE(m2 >= m2);

   REQUIRE(m1 == m1);
   REQUIRE(m2 == m1);
   REQUIRE(m1 == m2);
   REQUIRE(m2 == m2);

   REQUIRE(m1 != n1);
   REQUIRE(m2 != n1);
   REQUIRE(m1 != n2);
   REQUIRE(m2 != n2);
}

TEST_CASE("Verify integer power.", "[units]")
{
   length const x1 = 2 * km;
   dyndim const x2 = 2 * km;
   REQUIRE(x1.pow<2>() == 4 * km * km);
   REQUIRE(x1.pow<-2>() * (km * km) == Approx(0.25));
   REQUIRE(pow<2>(x1) == 4 * km * km);
   REQUIRE(pow<-2>(x1) * (km * km) == Approx(0.25));
   REQUIRE(x2.pow<2>() == 4 * km * km);
   REQUIRE((x2.pow<-2>() * (km * km)).number() == Approx(0.25));
   REQUIRE(x2.pow(2) == 4 * km * km);
   REQUIRE((x2.pow(-2) * (km * km)).number() == Approx(0.25));
   REQUIRE(pow<2>(x2) == 4 * km * km);
   REQUIRE((pow<-2>(x2) * (km * km)).number() == Approx(0.25));
   REQUIRE(pow(x2, 2) == 4 * km * km);
   REQUIRE((pow(x2, -2) * (km * km)).number() == Approx(0.25));
}

TEST_CASE("Verify integer root.", "[units]")
{
   volume const v1 = 8 * km * km * km;
   dyndim const v2 = 8 * km * km * km;
   REQUIRE(v1.root<3>() / km == Approx(2.0));
   REQUIRE(root<3>(v1) / km == Approx(2.0));
   REQUIRE(sqrt(v1 * 2 * km) / (km * km) == Approx(4.0));
   REQUIRE((v2.root<3>() / km).number() == Approx(2.0));
   REQUIRE((v2.root(3) / km).number() == Approx(2.0));
   REQUIRE((root<3>(v2) / km).number() == Approx(2.0));
   REQUIRE((root(v2, 3) / km).number() == Approx(2.0));
   REQUIRE((sqrt(v2 * 2 * km) / (km * km)).number() == Approx(4.0));
}

TEST_CASE("Verify absolute value.", "[units]")
{
   length const x1 = -1 * m;
   dyndim const x2 = -1 * m;
   REQUIRE(fabs(x1) == 1 * m);
   REQUIRE(fabs(x2) == 1 * m);
}

