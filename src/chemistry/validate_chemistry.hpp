/** @file validate_chemistry.hpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Thursday Jun 17, 2021 15:00:15 PDT
 * @bug No known bugs.
 */

#ifndef VALIDATE_CHEMISTRY_HPP
#define VALIDATE_CHEMISTRY_HPP

inline void validate_chemistry(Real const c[], Real const c0[], Real const c1[])
{
  for (int n = 0; n <= NVAPOR; ++n) {
    if (!(c[n] >= 0.)) {
      for (int t = 0; t < NHYDRO; ++t)
        std::cout << c[t] << " ";
      std::cout << std::endl;
      for (int t = 0; t < NHYDRO; ++t)
        std::cout << c0[t] << " ";
      std::cout << std::endl;
      for (int t = 0; t < NHYDRO; ++t)
        std::cout << c1[t] << " ";
      std::cout << std::endl << std::endl;
    }
    assert(c[n] >= 0.);
  }

  assert(c[IPR] >= 0.);
}

#endif /* end of include guard VALIDATE_CHEMISTRY_HPP */

