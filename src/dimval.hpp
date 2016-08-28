
#ifndef NUMERIC_DIMVAL_HPP
#define NUMERIC_DIMVAL_HPP

#include <cmath>   // for M_PI
#include <utility> // for std::forward()

namespace numeric
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
};

typedef dimval<1, 0, 0, 0, 0> time;
typedef dimval<0, 1, 0, 0, 0> length;
typedef dimval<0, 0, 1, 0, 0> mass;
typedef dimval<0, 0, 0, 1, 0> charge;
typedef dimval<0, 0, 0, 0, 1> temperature;
typedef dimval<-1, 1, 0, 0, 0> speed;
typedef dimval<-2, 1, 0, 0, 0> acceleration;
typedef dimval<-1, 1, 1, 0, 0> momentum;
typedef dimval<-2, 1, 1, 0, 0> force;
typedef dimval<-2, 2, 1, 0, 0> energy;
typedef dimval<-3, 2, 1, 0, 0> power;
typedef dimval<-1, 0, 0, 1, 0> current;

struct secs : public time {
   secs(double n = 0.0) { v_ = n; }
};

struct millisecs : public time {
   millisecs(double n = 0.0) { v_ = 0.001 * n; }
};

struct microsecs : public time {
   microsecs(double n = 0.0) { v_ = 1.0E-06 * n; }
};

struct nanosecs : public time {
   nanosecs(double n = 0.0) { v_ = 1.0E-09 * n; }
};

struct meters : public length {
   meters(double n = 0.0) { v_ = n; }
};

struct kilometers : public length {
   kilometers(double n = 0.0) { v_ = 1000.0 * n; }
};

struct millimeters : public length {
   millimeters(double n = 0.0) { v_ = 0.001 * n; }
};

struct micrometers : public length {
   micrometers(double n = 0.0) { v_ = 1.0E-06 * n; }
};

struct nanometers : public length {
   nanometers(double n = 0.0) { v_ = 1.0E-09 * n; }
};

struct angstroms : public length {
   angstroms(double n = 0.0) { v_ = 1.0E-10 * n; }
};

struct coulombs : public charge {
   coulombs(double n = 0.0) { v_ = n; }
};

struct kelvins : public temperature {
   kelvins(double n = 0.0) { v_ = n; }
};

struct newtons : public force {
   newtons(double n = 0.0) { v_ = n; }
};

struct dynes : public force {
   dynes(double n = 0.0) { v_ = 1.0E-05 * n; }
};

struct joules : public energy {
   joules(double n = 0.0) { v_ = n; }
};

struct ergs : public energy {
   ergs(double n = 0.0) { v_ = 1.0E-07 * n; }
};

struct watts : public power {
   watts(double n = 0.0) { v_ = n; }
};

struct amperes : public current {
   amperes(double n = 0.0) { v_ = n; }
};

extern time const sec;
extern time const millisec;
extern time const microsec;
extern time const nanosec;
extern length const meter;
extern length const kilometer;
extern length const millimeter;
extern length const micrometer;
extern length const nanometer;
extern length const angstrom;
extern charge const coulomb;
extern temperature const kelvin;
extern force const newton;
extern force const dyne;
extern energy const joule;
extern energy const erg;
extern power const watt;
extern current const ampere;
}

#endif // ndef NUMERIC_DIMVAL_HPP

