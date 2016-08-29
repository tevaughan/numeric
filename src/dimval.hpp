
#ifndef NUMERIC_DIMVAL_HPP
#define NUMERIC_DIMVAL_HPP

#include <cmath>    // for M_PI
#include <iostream> // for ostream
#include <utility>  // for forward()

namespace num
{
double constexpr deg = M_PI / 180.0;
double constexpr arcmin = deg / 60.0;
double constexpr arcsec = arcmin / 60.0;

inline double degs(double n) { return n * deg; }
inline double arcmins(double n) { return n * arcmin; }
inline double arcsecs(double n) { return n * arcsec; }

/// Model of a dimensioned value.
/// \tparam TI  Exponent of time.
/// \tparam D   Exponent of distance.
/// \tparam M   Exponent of mass.
/// \tparam C   Exponent of charge.
/// \tparam TE  Exponent of temperature.
template <int TI, int D, int M, int C, int TE> class dimval
{
   /// Used by friend integral().
   operator double() const { return v_; }

protected:
   double v_; ///< Value in MKS.

   dimval(double vv) : v_(vv) {}

public:
   dimval() : v_(0.0) {}

   template <int OTI, int OD, int OM, int OC, int OTE> friend class dimval;

   template <typename F, typename A>
   friend auto integral(F f, A a, A b, double t, unsigned n)
         -> decltype(std::forward<F>(f)(a) * a);

   /// Add dimensioned values.
   dimval operator+(dimval dv) const { return v_ + dv; }

   /// Subtract dimensioned values.
   dimval operator-(dimval dv) const { return v_ - dv; }

   /// Additive assignment.
   dimval &operator+=(dimval dv) { return v_ += dv; }

   /// Subtractive assignment.
   dimval &operator-=(dimval dv) { return v_ -= dv; }

   /// Multiply dimensioned values.
   template <int OTI, int OD, int OM, int OC, int OTE>
   dimval<TI + OTI, D + OD, M + OM, C + OC, TE + OTE>
   operator*(dimval<OTI, OD, OC, OM, OTE> dv) const
   {
      return v_ * dv;
   }

   /// Multiplication by number on right side.
   dimval operator*(double s) const { return v_ * s; }

   /// Multiplication of dimval by number on left side.
   friend dimval operator*(double s, dimval dv) { return s * dv.v_; }

   /// Division resulting in double.
   double operator/(dimval dv) const { return v_ / dv; }

   /// Divide dimensioned values.
   template <int OTI, int OD, int OM, int OC, int OTE>
   dimval<TI - OTI, D - OD, M - OM, C - OC, TE - OTE>
   operator/(dimval<OTI, OD, OM, OC, OTE> dv) const
   {
      return v_ / dv;
   }

   /// Division by number.
   dimval operator/(double s) const { return v_ / s; }

   /// Multiplicative assignment against double.
   dimval &operator*=(double s) { return v_ *= s; }

   /// Multiplicative division against double.
   dimval &operator/=(double s) { return v_ /= s; }

   bool operator<(dimval dv) const { return v_ < dv; }
   bool operator>(dimval dv) const { return v_ > dv; }
   bool operator<=(dimval dv) const { return v_ <= dv; }
   bool operator>=(dimval dv) const { return v_ >= dv; }
   bool operator==(dimval dv) const { return v_ == dv; }
   bool operator!=(dimval dv) const { return v_ != dv; }

   /// Exponentiation.
   template <int E> dimval<TI * E, D * E, M * E, C * E, TE * E> pow() const
   {
      return pow(v_, E);
   }

   /// Root.
   template <unsigned E>
   dimval<TI / E, D / E, M / E, C / E, TE / E> pow() const
   {
      static_assert(E, "zeroth root");
      static_assert(D / E * E == D, "illegal root along distance");
      static_assert(M / E * E == M, "illegal root along mass");
      static_assert(C / E * E == C, "illegal root along charge");
      static_assert(TI / E * E == TI, "illegal root along time");
      static_assert(TE / E * E == TE, "illegal root along temperature");
      double constexpr e = 1.0 / E;
      return pow(v_, e);
   }

   friend dimval fabs(dimval dv) { return fabs(dv.v_); }

   friend std::ostream &operator<<(std::ostream &os, dimval dv)
   {
      os << "[" << dv.v_;
      if (M == 1) {
         os << " kg";
      } else if(M != 0) {
         os << " kg^" << M;
      }
      if (D == 1) {
         os << " m";
      } else if(D != 0) {
         os << " m^" << D;
      }
      if (TI == 1) {
         os << " s";
      } else if(TI != 0) {
         os << " s^" << TI;
      }
      if (C == 1) {
         os << " C";
      } else if(C != 0) {
         os << " C^" << C;
      }
      if (TE == 1) {
         os << " K";
      } else if(TE != 0) {
         os << " K^" << TE;
      }
      return os << "]";
   }
};

/// Specialization for dimensionless value.
template <> class dimval<0, 0, 0, 0, 0>
{
public:
   double v;
   dimval(double vv = 0.0) : v(vv) {}
   template <int OTI, int OD, int OM, int OC, int OTE> friend class dimval;
   operator double const &() const { return v; }
   operator double &() { return v; }

   /// Add dimensioned values.
   double operator+(dimval dv) const { return v + dv; }

   /// Subtract dimensioned values.
   double operator-(dimval dv) const { return v - dv; }

   /// Additive assignment.
   double &operator+=(dimval dv) { return v += dv; }

   /// Subtractive assignment.
   double &operator-=(dimval dv) { return v -= dv; }

   /// Multiply dimensioned values.
   template <int OTI, int OD, int OM, int OC, int OTE>
   dimval<OTI, OD, OM, OC, OTE>
   operator*(dimval<OTI, OD, OC, OM, OTE> dv) const
   {
      return v * dv;
   }

   /// Multiplication by number on right side.
   double operator*(double s) const { return v * s; }

   /// Multiplication of dimval by number on left side.
   friend double operator*(double s, dimval dv) { return s * dv.v; }

   /// Division resulting in double.
   double operator/(dimval dv) const { return v / dv; }

   /// Divide dimensioned values.
   template <int OTI, int OD, int OM, int OC, int OTE>
   dimval<-OTI, -OD, -OM, -OC, -OTE>
   operator/(dimval<OTI, OD, OM, OC, OTE> dv) const
   {
      return v / dv;
   }

   /// Division by number.
   double operator/(double s) const { return v / s; }

   /// Multiplicative assignment against double.
   double &operator*=(double s) { return v *= s; }

   /// Multiplicative division against double.
   double &operator/=(double s) { return v /= s; }

   bool operator<(dimval dv) const { return v < dv; }
   bool operator>(dimval dv) const { return v > dv; }
   bool operator<=(dimval dv) const { return v <= dv; }
   bool operator>=(dimval dv) const { return v >= dv; }
   bool operator==(dimval dv) const { return v == dv; }
   bool operator!=(dimval dv) const { return v != dv; }

   bool operator<(double dv) const { return v < dv; }
   bool operator>(double dv) const { return v > dv; }
   bool operator<=(double dv) const { return v <= dv; }
   bool operator>=(double dv) const { return v >= dv; }
   bool operator==(double dv) const { return v == dv; }
   bool operator!=(double dv) const { return v != dv; }

   /// Exponentiation.
   template <int E> double pow() const { return std::pow(v, E); }

   /// Root.
   template <unsigned E> double pow() const
   {
      double constexpr e = 1.0 / E;
      return std::pow(v, e);
   }

   friend double fabs(dimval dv) { return fabs(dv.v); }

   friend std::ostream &operator<<(std::ostream &os, dimval dv)
   {
      return os << dv.v;
   }
};
}

#endif // ndef NUMERIC_DIMVAL_HPP

