#ifndef CELL_VARIABLES_HPP
#define CELL_VARIABLES_HPP

// C/C++ header
#include <vector>

// Athena++ header
#include "../athena.hpp"

//! \class CellVariables
//  \brief a collection of all physical data in a computational cell
class CellVariables {
public:
  //! data pointers
  //! hydro data
  Real *q;
  //! scalar data
  Real *s;
  //! particle data
  Real *c;

  // functions
  CellVariables(int num) :
    data_(num)
  {
    q = data_.data();
    s = data_.data() + NHYDRO;
    c = data_.data() + NHYDRO + NSCALARS;
  }

private:
  // data holder
  std::vector<Real> data_;
};

#endif
