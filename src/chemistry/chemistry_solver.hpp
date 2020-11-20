#ifndef CHEMISTRY_SOLVER_HPP
#define CHEMISTRY_SOLVER_HPP

//#include "../math/linalg.h"
// Athena++ headers
#include "../athena.hpp"

// Eigen header files
#include "../math/eigen335/Eigen/Core"
#include "../math/eigen335/Eigen/Dense"
#include "../math/eigen335/Eigen/LU"

class ChemistrySolver {
public:
  // needed for Eigen small matrix
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  ChemistrySolver() {}

  // solve independent chemistry: A*dq = b
  // the second argument Real *dq is both input (b) and output (dq)
  void SolveIndependent(Real **A, Real *dq) {
    for (int i = 1; i <= NVAPOR; ++i) {
      // assemble matrix for each species
      #pragma ivdep
      for (int j = 0; j < NPHASE; ++j) {
        for (int k = 0; k < NPHASE; ++k)
          a_(j,k) = A[j*NVAPOR+i][k*NVAPOR+i];
        sol_(j) = dq[j*NVAPOR+i];
      }

      //ludcmp<NPHASE>(a_, indx_, vv_);
      //lubksb<NPHASE>(a_, indx_, sol_);

      if (NPHASE <= 4) sol_ = a_.inverse()*sol_;
      else sol_ = a_.partialPivLu().solve(sol_);

      for (int j = 0; j < NPHASE; ++j)
        dq[j*NVAPOR+i] = sol_(j);
    }
    dq[0] = 0.;
  }

private:
  // scratch arrays
  Eigen::Matrix<Real, NPHASE, NPHASE> a_;
  Eigen::Matrix<Real, NPHASE, 1> sol_;

  //Real a_[NPHASE][NPHASE];
  //Real sol_[NPHASE], vv_[NPHASE];
  //int indx_[NPHASE];
};

#endif
