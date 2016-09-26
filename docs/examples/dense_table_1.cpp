
#include <iostream>
#include "dense-table.hpp"

using namespace num;
using namespace std;

class F
{
   double v_;

public:
   F(double vv = 0.0) : v_(vv) {}
   double operator()(double) const { return v_; }
};

int main()
{
   std::vector<F> vf;
   vf.push_back(F(1.0));
   vf.push_back(F(2.0));
   dense_table<double, F> t(3.0, 1.1, vf);
   for (double x = 0.0; x < 5.0; x += 0.05) {
      cout << x << " " << t(x) << endl;
   }
   return 0;
}

