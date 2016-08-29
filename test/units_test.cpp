
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
   temperature const t2 = 2.5 * μK;
   ostringstream oss3, oss4;
   oss3 << t1;
   oss4 << t2;
   REQUIRE(oss3.str() == "[2.5 μK]");
   REQUIRE(oss4.str() == "[2.5e-06 K]");
}

