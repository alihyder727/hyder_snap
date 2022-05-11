// C/C++ header
#include "grid_data.hpp"

GridData::GridData() {
  q = data_;
  s = data_ + NHYDRO;
  c = data_ + NHYDRO + NSCALARS;
}
