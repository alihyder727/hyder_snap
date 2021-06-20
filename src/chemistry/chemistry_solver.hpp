#ifndef CHEMISTRY_SOLVER_HPP
#define CHEMISTRY_SOLVER_HPP

//#include "../math/linalg.h"
// Athena++ headers
#include "../athena.hpp"

// Eigen header files
#include "../math/eigen335/Eigen/Core"
#include "../math/eigen335/Eigen/Dense"
#include "../math/eigen335/Eigen/LU"

template<int N>
class ChemistrySolver {
public:
// needed for Eigen small matrix
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

// data
  enum {Size = N};

// functions
  ChemistrySolver() {
    I_.setIdentity();
  }

  template<typename T1, typename T2>
  T1 BDF1(T1 const& Rate, T2 const& Jac, Real dt) {
    A_ = I_/dt - Jac;
    if (N <= 4) return A_.inverse()*Rate;
    else return A_.partialPivLu().solve(Rate);
  }

  template<typename T1, typename T2>
  T1 TRBDF2(T1 const& Rate, T2 const& Jac, Real dt) {
    int gamma = 2. - sqrt(2.);
    A_ = I_ - gamma/2.*dt*Jac;
    B_ = dt*(1. - gamma/2.)*Rate + dt*gamma/2.*A_*Rate;
    A_ = A_*A_;
    if (N <= 4) return A_.inverse()*B_;
    else return A_.partialPivLu().solve(B_);
  }

private:
  // scratch array
  Eigen::Matrix<Real,N,N> A_, I_;
  Eigen::Matrix<Real,N,1> B_;
};

#endif
