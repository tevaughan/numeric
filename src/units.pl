#!/usr/bin/env perl

# Copyright 2016  Thomas E. Vaughan
# Distributable under the terms of the GNU LGPL, Version 3 or later.

open(INPUT, '<', "units.txt");
open(OUT_H, '>', "units.hpp");
open(OUT_C, '>', "units.cpp");

print OUT_H <<'EOF1';

// This file was generated by the 'units' executable.

/// \file   units.hpp
/// \brief  Automatically generated header file for units.

#ifndef NUMERIC_UNITS_HPP
#define NUMERIC_UNITS_HPP

#include "dimensions.hpp"

namespace num {
EOF1

print OUT_C << 'EOF2';

// This file was generated by the 'units' executable.

/// \file   units.cpp
/// \brief  Automatically generated implementation file for units.

#include "units.hpp"

using namespace num;

EOF2

while(<INPUT>) {
   s/#.*$//;         # Remove comment.
   next if /^\s*$/;  # Skip blank line.
   my($pnm, $snm, $dim, $cnv) = split;
   print OUT_H <<"EOF3";

/// Machine-generated structure providing public constructor, which takes a
/// double-precision number of $pnm.
struct $pnm : public $dim {
   /// Construct a dimensioned quantity from a double-precision number of
   /// $pnm.
   $pnm(double v) : $dim(v * $cnv) {}
   /// Allow default copying.
   $pnm($pnm const& dv) = default;
   /// Write representation to output stream.
   friend std::ostream& operator<<(std::ostream& os, $pnm const& dv)
   {
      /// Factor that converts to the number of $pnm from the number of the
      /// equivalent MKS unit.
      static double constexpr c = 1.0 / $cnv;
      return os << "[" << c * dv.v_ << " $snm]";
   }
};

namespace u {
/// Declaration of symbol $snm representing a unit of $dim.
extern $pnm const $snm;
}
EOF3
   print OUT_C "$pnm const num::u::$snm(1.0);\n";
}

print OUT_H <<'EOF4';
}

#endif // ndef NUMERIC_UNITS_HPP

EOF4

print OUT_C "\n";

