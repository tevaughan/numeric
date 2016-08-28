
#include "catch.hpp"
#include "integral.hpp"

using namespace num;
using namespace std;

TEST_CASE("Verify numeric integration of linear function.", "[integral]")
{
   auto lambda = [](double x) { return x; };
   REQUIRE(integral(lambda, 0.0, 1.0) == Approx(0.5));
}

