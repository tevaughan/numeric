
#include <sstream> // for ostringstream

#include "catch.hpp"
#include "integral.hpp"
#include "units.hpp"

using namespace num;
using namespace std;

TEST_CASE("Verify integration of lambda.", "[integral]")
{
   function<double(double)> linear = [](double x) { return x; };
   REQUIRE(integral(linear, 0, 1) == Approx(0.5));
}

area square(length x) { return x * x; }

TEST_CASE("Verify integration of (dimval) global function.", "[integral]")
{
   // Verify that integral of x^2 from  0 cm  to  1 cm  is  1/3 cm^3.
   volume const i = integral(square, 0 * cm, 1 * cm);
   ostringstream oss;
   oss << i;
   REQUIRE(oss.str() == "[3.33333e-07 m^3]");
   REQUIRE(i / cm.pow<3>() == Approx(1.0 / 3.0));
}

struct foo {
   double sin(double x) { return std::sin(x); }
};

TEST_CASE("Verify integration of member function.", "[integral]")
{
   foo f;
   function<double(double)> s = [&f](double x) { return f.sin(x); };
   double const i = integral(s, 0, M_PI);
}

