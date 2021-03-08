// C/C++ headers
#include <cmath>

// Athena++ header
#include "../athena.hpp"

int IsPositive(Real v) {
  return v > 0 ? 1 : 0;
}

int IsNumber(Real v) {
  return !std::isnan(v);
}
