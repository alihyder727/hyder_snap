// Athena++ header
#include "../athena.hpp"

//! \class GridData
//  \brief a collection of all physical data in a computational cell
class GridData {
public:
  // data
  Real *q;
  Real *s;
  Real *c;

  // functions
  GridData();

private:
  Real data_[NGRIDMAX];
};
