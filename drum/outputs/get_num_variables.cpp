/** \file get_num_variables.cpp
 * \brief provides function definition, getNumVariables
 *
 * \author Cheng Li (chengcli@umich.edu)
 * \date Monday Apr 25, 2022 11:58:37 EDT
 */

#include "outputs.hpp"

int getNumVariables(std::string grid, AthenaArray<Real> const& data)
{
  int nvar;
  if (grid == "--C" || grid == "--F") {
    nvar = data.GetDim2();
  } else if (grid == "---") {
    nvar = data.GetDim1();
  } else {
    nvar = data.GetDim4();
  }

  return nvar;
}
