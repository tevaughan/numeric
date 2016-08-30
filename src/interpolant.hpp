
#ifndef NUMERIC_INTERPOLANT_HPP
#define NUMERIC_INTERPOLANT_HPP

#include <algorithm> // for sort()
#include <fstream>   // for ifstream
#include <memory>    // for unique_ptr
#include <sstream>   // for istringstream
#include <string>    // for string
#include <utility>   // for pair
#include <vector>    // for vector

namespace num
{
   /// Point used as one of the constraints of a linear interpolant.
   /// \tparam I  Type of independent variable (x value).
   /// \tparam D  Type of dependent variable (y value).
   template <typename I, typename D>
   using ipoint = std::pair<I, D>;

   /// List of points used to constrain a linear interpolant.
   ///
   /// An instance of ilist can be used to initialize an instance of
   /// interpolant.
   ///
   /// \tparam I  Type of independent variable (x value).
   /// \tparam D  Type of dependent variable (y value).
   template <typename I, typename D>
   using ilist = std::vector<ipoint<I, D>>;

   /// Linear interpolant.
   ///
   /// This original C++ implementation is due entirely to Thomas E. Vaughan,
   /// who wrote it in 2016.
   ///
   /// Class interpolant is designed with the idea that either of its template
   /// type parameters may be a double or else a plural-name descendant (like
   /// meters) of dimval.  In that way, the corresponding value can always be
   /// constructed from an initializer of type double.
   ///
   /// \tparam I  Type of independent variable (x value).
   /// \tparam D  Type of dependent variable (y value).
   template <typename I, typename D>
   class interpolant
   {
      using list = ilist<I, D>;
      using point = ipoint<I, D>;

      list data_;

      static bool cmpr_pts(point const &p1, point const &p2)
      {
         return p1.first < p2.first;
      }

      void sort() { std::sort(data_.begin(), data_.end(), cmpr_pts); }

      static std::unique_ptr<std::ifstream> ifstr(std::string fname)
      {
         using namespace std;
         unique_ptr<ifstream> p(new ifstream(fname));
         if (!(*p)) {
            throw "failed to open '" + fname + "'";
         }
         return p;
      }

      static point get_point(std::string line, I x_unit, D y_unit)
      {
         using namespace std;
         istringstream iss(line);
         double x, y;
         if (!(iss >> x)) {
            throw "error reading independent variable";
         }
         if (!(iss >> y)) {
            throw "error reading dependent variable";
         }
         return point(x * x_unit, y * y_unit);
      }

   public:
      /// Initialize ifrom an ilist.
      /// \param d  Data stored in ilist.
      interpolant(list d = list()) : data_(std::move(d)) { sort(); }

      /// Initialize from the first two columns of a space-delimited ASCII
      /// file.
      /// \param fname  Name of input file.
      interpolant(std::string fname, I x_unit = 1.0, D y_unit = 1.0)
      {
         using namespace std;
         unique_ptr<ifstream> is = ifstr(fname);
         string line;
         while (getline(*is, line)) {
            size_t const p = line.find_first_not_of(" \f\n\r\t");
            if (p == string::npos || line[p] == '#') {
               continue;
            }
            data_.push_back(get_point(line, x_unit, y_unit));
         }
         sort();
      }

      /// Interpolate.
      /// \param x  Value of independent variable.
      /// \return   Interpolated value of dependent variable.
      D operator()(I x) const
      {
         if (data_.size() == 0) {
            throw "interpolating on empty interpolant";
         }
         if (x <= data_.begin()->first) {
            return data_.begin()->second;
         }
         if (x >= data_.rbegin()->first) {
            return data_.rbegin()->second;
         }
         using namespace std;
         point const xp{x, data_[0].second}; // Dummy y-coord used for search.
         auto const j = upper_bound(data_.begin(), data_.end(), xp, cmpr_pts);
         auto const i = j - 1;
         I const &xi = i->first;
         I const &xj = j->first;
         D const &yi = i->second;
         D const &yj = j->second;
         return yi + (yj - yi) * ((x - xi) / (xj - xi));
      }
   };

   using interpolantd = interpolant<double, double>;
   using ilistd = ilist<double, double>;
}

#endif // ndef NUMERIC_INTERPOLANT_HPP

