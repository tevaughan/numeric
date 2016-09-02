
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
   template <int TI, int D, int M, int C, int TE>
   class dimval;

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
   /// type parameters may be a double or else a dimval.
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
      ///
      /// \param fname   Name of input file.
      ///
      /// \param x_unit  Factor multiplied against each number in file's first
      ///                space-delimited column. Product is stored as x value.
      ///
      /// \param y_unit  Factor multiplied against each number in file's second
      ///                space-delimited column. Product is stored as y value.
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

      /// Multiply two interpolants together such that there is a point in the
      /// product interpolant at every unique x coordinate across both
      /// interpolants.  Each new y coordinate is the product of a y coordinate
      /// in one interpolant and the corresponding interpolated y coordinate in
      /// the other.
      ///
      /// \tparam OD    Type of y coordinate in other interpolant.
      /// \param  oint  Reference to other interpolant.
      /// \return       Product interpolant.
      template <typename OD>
      auto operator*(interpolant<I, OD> const &oint) const
            -> interpolant<I, decltype(data_[0].second *oint(data_[0].first))>
      {
         using PD = decltype(data_[0].second *oint(data_[0].first));
         ilist<I, PD> pl;
         ilist<I, D> const& l1 = data_;
         ilist<I, OD> const& l2 = oint.data_;
         unsigned const s1 = l1.size();
         unsigned const s2 = l2.size();
         unsigned n1 = 0;
         unsigned n2 = 0;
         I ind;
         if (s1 == 0 || s2 == 0) {
            return pl;
         }
         while (n1 < s1 || n2 < s2) {
            if (n1 < s1 && n2 < s2) {
               I const i1 = l1[n1].first;
               I const i2 = l2[n2].first;
               D const d1 = l1[n1].second;
               OD const d2 = l2[n2].second;
               if (i1 < i2) {
                  pl.push_back({i1, d1 * oint(i1)});
                  ++n1;
               } else {
                  pl.push_back({i2, (*this)(i2) * d2});
                  ++n2;
               }
            } else if (n1 < s1) {
               I const i1 = l1[n1].first;
               D const d1 = l1[n1].second;
               pl.push_back({i1, d1 * oint(i1)});
               ++n1;
            } else {
               I const i2 = l2[n2].first;
               OD const d2 = l2[n2].second;
               pl.push_back({i2, (*this)(i2)*d2});
               ++n2;
            }
         }
         return pl;
      }

      /// Multiply every y-value of interpolant by a scale factor on the
      /// right, and return the resultant interpolant.
      ///
      /// \param  s  Scale factor for y-values.
      /// \return    Scaled interpolant.
      auto operator*(double s) const
            -> interpolant<I, decltype(data_[0].second *s)>
      {
         using ND = decltype(data_[0].second *s);
         ilist<I, ND> nl(data_.size());
         for (unsigned i = 0; i < data_.size(); ++i) {
            nl[i] = ipoint<I, ND>(data_[i].first, data_[i].second * s);
         }
         return nl;
      }

      /// Multiply every y-value of interpolant by a scale factor on the left,
      /// and return the resultant interpolant.  Assume that multiplication of
      /// y-values is commutative.
      ///
      /// \param  s  Scale factor for y-values.
      /// \param  i  Interpolant to scale.
      /// \return    Scaled interpolant.
      friend auto operator*(double s, interpolant const &i)
            -> interpolant<I, decltype(data_[0].second *s)>
      {
         return i * s;
      }

      /// Multiply every y-value of interpolant by a dimval scale factor on the
      /// right, and return the resultant interpolant.
      ///
      /// \param  s  Scale factor for y-values.
      /// \return    Scaled interpolant.
      template <int TI, int ID, int M, int C, int TE>
      auto operator*(dimval<TI, ID, M, C, TE> s) const
            -> interpolant<I, decltype(data_[0].second *s)>
      {
         using ND = decltype(data_[0].second *s);
         ilist<I, ND> nl(data_.size());
         for (unsigned i = 0; i < data_.size(); ++i) {
            nl[i] = ipoint<I, ND>(data_[i].first, data_[i].second * s);
         }
         return nl;
      }

      /// Multiply every y-value of interpolant by a dimval scale factor on the
      /// left, and return the resultant interpolant.  Assume that
      /// multiplication of y-values is commutative.
      ///
      /// \param  s  Scale factor for y-values.
      /// \param  i  Interpolant to scale.
      /// \return    Scaled interpolant.
      template <int TI, int ID, int M, int C, int TE>
      friend auto operator*(dimval<TI, ID, M, C, TE> s, interpolant const &i)
            -> interpolant<I, decltype(data_[0].second *s)>
      {
         return i * s;
      }

      /// Divide every y-value of interpolant by a scale factor, and return the
      /// resultant interpolant.
      ///
      /// \tparam Y  Type of scale factor for y-values.
      /// \param  s  Scale factor for y-values.
      /// \return    Scaled interpolant.
      template <typename Y>
      auto operator/(Y s) const
            -> interpolant<I, decltype(data_[0].second / s)>
      {
         using ND = decltype(data_[0].second / s);
         ilist<I, ND> nl(data_.size());
         for (unsigned i = 0; i < data_.size(); ++i) {
            nl[i] = ipoint<I, ND>(data_[i].first, data_[i].second / s);
         }
         return nl;
      }
   };

   using interpolantd = interpolant<double, double>;
   using ilistd = ilist<double, double>;
}

#endif // ndef NUMERIC_INTERPOLANT_HPP

