
#include <sstream> // for ostringstream

#include "catch.hpp"
#include "integral.hpp"
#include "units.hpp"

using namespace num;
using namespace std;

TEST_CASE("Verify integration of lambda.", "[integral]")
{
   function<double(double)> linear = [](double x) { return x; };
   // Verify that integral of x of 0 to 1 is 1/2.
   REQUIRE(integral(linear, 0, 1) == Approx(0.5));
}

area square(length x) { return x * x; }

TEST_CASE("Verify integration of (dimval) global function.", "[integral]")
{
   volume const i = integral(square, 0 * cm, 1 * cm);
   ostringstream oss;
   oss << i;
   // Verify that integral of x^2 from  0 cm  to  1 cm  is  1/3 cm^3.
   REQUIRE(oss.str() == "[3.33333e-07 m^3]");
   REQUIRE(i / cm.pow<3>() == Approx(1.0 / 3.0));
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
   double const i = integral(s, 0, M_PI);
   // Verify that integral of sin(x) from 0 to pi is 2.
   REQUIRE(i == Approx(2.0));
   f.freq = 2.0;
   double const j = integral(s, 0, M_PI);
   // Verify that integral of sin(2*x) from 0 to pi is 0.
   REQUIRE(j == Approx(0.0));
}

TEST_CASE("Trigger coverage of code requiring at least two samples.",
          "[integral]")
{
   volume const i = integral(square, 1 * cm, 2 * cm, 1.0E-06, 0);
   // Verify that integral of x^2 from 1 cm  to  2 cm  is  7/3 cm^3.
   REQUIRE(i / cm.pow<3>() == Approx(7.0 / 3.0));
}

