#ifndef FORWARD_BACKWARD_HPP
#define FORWARD_BACKWARD_HPP

// C/C++ headers
#include <vector>

// Eigen headers
#include "../../math/eigen335/Eigen/Core"
#include "../../math/eigen335/Eigen/Dense"

// Athena++ headers
#include "implicit_solver.hpp"

template<typename T1, typename T2>
void ImplicitSolver::ForwardSweep(
  std::vector<T1> &a, std::vector<T1> &b, std::vector<T1> &c, 
  std::vector<T2> &delta, std::vector<T2> &corr, Real dt,
  int k, int j, int il, int iu)
{
  T2 rhs;

  if (T2::RowsAtCompileTime == 3) {  // partial matrix
    rhs(0) = du_(IDN,k,j,il)/dt;
    rhs(1) = du_(IVX+mydir,k,j,il)/dt;
    rhs(2) = du_(IEN,k,j,il)/dt;
    rhs -= corr[il];
  } else {  // full matrix
    rhs(0) = du_(IDN,k,j,il)/dt;
    rhs(1) = du_(IVX+mydir,k,j,il)/dt;
    rhs(2) = du_(IVX+(IVY-IVX+mydir)%3,k,j,il)/dt;
    rhs(3) = du_(IVX+(IVZ-IVX+mydir)%3,k,j,il)/dt;
    rhs(4) = du_(IEN,k,j,il)/dt;
  }

  if (has_bot_neighbor) {
    RecvBotBuffer(a[il-1], delta[il-1], k, j, bblock);
    a[il] = (a[il] - b[il]*a[il-1]).inverse().eval();
    delta[il] = a[il]*(rhs - b[il]*delta[il-1]);
    a[il] *= c[il];
  } else {
    a[il] = a[il].inverse().eval();
    delta[il] = a[il]*rhs;
    a[il] *= c[il];
  }

  for (int i = il+1; i <= iu; ++i) {
    if (T2::RowsAtCompileTime == 3) {  // partial matrix
      rhs(0) = du_(IDN,k,j,i)/dt;
      rhs(1) = du_(IVX+mydir,k,j,i)/dt;
      rhs(2) = du_(IEN,k,j,i)/dt;
      rhs -= corr[i];
    } else {
      rhs(0) = du_(IDN,k,j,i)/dt;
      rhs(1) = du_(IVX+mydir,k,j,i)/dt;
      rhs(2) = du_(IVX+(IVY-IVX+mydir)%3,k,j,i)/dt;
      rhs(3) = du_(IVX+(IVZ-IVX+mydir)%3,k,j,i)/dt;
      rhs(4) = du_(IEN,k,j,i)/dt;
    }

    a[i] = (a[i] - b[i]*a[i-1]).inverse().eval();
    delta[i] = a[i]*(rhs - b[i]*delta[i-1]);
    a[i] *= c[i];
  }

  SaveCoefficients(a, delta, k, j, il, iu);

  if (has_top_neighbor)
    SendTopBuffer(a[iu], delta[iu], k, j, tblock);
}

template<typename T1, typename T2>
void ImplicitSolver::BackwardSubstitution(
  std::vector<T1> &a, 
  std::vector<T2> &delta, 
  int kl, int ku, int jl, int ju, int il, int iu)
{
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j) {
      LoadCoefficients(a, delta, k, j, il, iu);
      if (has_top_neighbor) {
        RecvTopBuffer(delta[iu+1], k, j, tblock);
        delta[iu] -= a[iu]*delta[iu+1];
      }

      // update solutions, i=iu
      for (int i = iu-1; i >= il; --i)
        delta[i] -= a[i]*delta[i+1];

      // 7. update conserved variables, i = iu
      for (int i = il; i <= iu; ++i) {
        if (T2::RowsAtCompileTime == 3) {  // partial matrix
          du_(IDN,k,j,i) = delta[i](0);
          du_(IVX+mydir,k,j,i) = delta[i](1);
          du_(IEN,k,j,i) = delta[i](2);
        } else { // full matrix
          du_(IDN,k,j,i) = delta[i](0);
          du_(IVX+mydir,k,j,i) = delta[i](1);
          du_(IVX+(IVY-IVX+mydir)%3,k,j,i) = delta[i](2);
          du_(IVX+(IVZ-IVX+mydir)%3,k,j,i) = delta[i](3);
          du_(IEN,k,j,i) = delta[i](4);
        }
      }

      if (has_bot_neighbor)
        SendBotBuffer(delta[il], k, j, bblock);
    }
}

#endif
