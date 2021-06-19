/** @file validate_chemistry.hpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Thursday Jun 17, 2021 15:00:15 PDT
 * @bug No known bugs.
 */

#ifndef VALIDATE_CHEMISTRY_HPP
#define VALIDATE_CHEMISTRY_HPP

template<typename T>
inline void validate_chemistry(Real const q[], AthenaArray<Real> const& c, 
  Real const q0[], Real const q1[], Real const q2[], int k, int j, int i, T rate)
{
  for (int n = 0; n <= NVAPOR; ++n) {
    if (!(q[n] >= 0.)) {
      std::cout << rate << std::endl;
      for (int t = 0; t < NHYDRO; ++t)
        std::cout << q[t] << " ";
      std::cout << std::endl;
      for (int t = 0; t < NHYDRO; ++t)
        std::cout << q0[t] << " ";
      std::cout << std::endl;
      for (int t = 0; t < NHYDRO; ++t)
        std::cout << q1[t] << " ";
      std::cout << std::endl;
      for (int t = 0; t < NHYDRO; ++t)
        std::cout << q2[t] << " ";
      std::cout << std::endl << std::endl;
    }
    assert(q[n] >= 0.);
  }

  assert(q[IPR] >= 0.);

  for (int t = 0; t < c.GetDim4(); ++t)
    assert(c(t,k,j,i) >= 0.);
}

#endif /* end of include guard VALIDATE_CHEMISTRY_HPP */

