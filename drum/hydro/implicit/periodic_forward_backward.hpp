#ifndef PERIODIC_FORWARD_BACKWARD_HPP
#define PERIODIC_FORWARD_BACKWARD_HPP

// C/C++ headers
#include <vector>

// Eigen headers
#include "../../math/eigen335/Eigen/Core"
#include "../../math/eigen335/Eigen/Dense"

// Athena++ headers
#include "communication.hpp"

template<typename T1, typename T2>
void ImplicitSolver::PeriodicForwardSweep(
  std::vector<T1> &diag, std::vector<T1> &diagL, std::vector<T1> &diagU,
  std::vector<T2> &delta, std::vector<T2> &corr, Real dt,
  int k, int j, int il, int iu)
{
  T1 corner, sum_beta_gamma;
  T2 rhs, sum_beta_zeta;
  std::vector<T1> alpha(diag.size()), beta(diag.size()), gamma(diag.size());
  std::vector<T1> ainv(diag.size());
  std::vector<T2> zeta(diag.size());

  // start periodic linear system solver, ref: El-Mikkawy (2005).
  // step 1: compute alpha, beta, gamma, delta from is to ie-1.
  //                                                 -1
  // alpha[1] = a[1], gamma[1] = U, beta[1] = L*alpha[1]

  // diag -> d
  // diagU -> a
  // diagL -> b
  // alpha -> c
  // gamma -> v
  // beta -> h
  // zeta -> k
  if (last_block)
    SendBuffer(diagU[iu], 0, 0, tblock);

  if (first_block)
    RecvBuffer(corner, 0, 0, bblock);
  else
    RecvBuffer(ainv[il-1], gamma[il-1], beta[il-1], zeta[il-1], 
      sum_beta_gamma, sum_beta_zeta, k, j, bblock);

  for (int i = il; i <= iu; ++i) {
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
    
    if ((i == il) && first_block) { // I'm the first row of the first block
      // c[1] = d[1]
      alpha[il] = diag[il]; 
      ainv[il] = alpha[il].inverse().eval();
      // v[1] = t (upper corner)
      gamma[il] = diagL[il];
      //             -1
      // h[1] = s*c[1]
      beta[il] = corner*ainv[il];
      // r[1] = k[1]
      zeta[il] = rhs;
      sum_beta_gamma.setZero();
      sum_beta_zeta.setZero();
      continue;
    }

    //                         -1
    // c[i] = d[i] - b[i]*c[i-1]*a[i-1]
    alpha[i] = diag[i] - diagL[i]*ainv[i-1]*diagU[i-1];
    ainv[i] = alpha[i].inverse().eval();
    //                   -1
    // v[i] = -b[i]*c[i-1]*v[i-1]
    gamma[i] = -diagL[i]*ainv[i-1]*gamma[i-1];
    //                   -1
    // h[i] = -a[i-1]*c[i]*h[i-1]
    beta[i] = -diagU[i-1]*ainv[i]*beta[i-1];
    //                       -1
    // r[i] = k[i] - b[i]*c[i]*r[i-1]
    zeta[i] = rhs - diagL[i]*ainv[i-1]*zeta[i-1];
  }


  if (last_block) {  // I'm the last block
    //                               -1
    // v[n-1] = a[n-1] - b[n-1]*c[n-2]*v[n-2];
    gamma[iu-1] = diagU[iu-1] - diagL[iu-1]*ainv[iu-2]*gamma[iu-2];
    //                                      -1
    // h[n-1] = (b[n] - h[n-2]*a[n-2])*c[n-1]
    beta[iu-1] = (diagL[iu] - beta[iu-2]*diagU[iu-2])*ainv[iu-1];
    for (int i = il; i < iu; ++i) {
      // h[i]*v[i]
      sum_beta_gamma += beta[i]*gamma[i];
      // h[i]*r[i]
      sum_beta_zeta += beta[i]*zeta[i];
    }
    //                 n-1
    // c[n] = d[n] - Sum {h[i]*v[i]}
    //                 i=1
    alpha[iu] = diag[iu] - sum_beta_gamma;
    ainv[iu] = alpha[iu].inverse().eval();
    //                 n-1
    // r[n] = k[n] - Sum {h[i]*r[i]}
    //                 i=1
    zeta[iu] = rhs - sum_beta_zeta;
  } else {
    for (int i = il; i <= iu; ++i) {
      // h[i]*v[i]
      sum_beta_gamma += beta[i]*gamma[i];
      // h[i]*r[i]
      sum_beta_zeta += beta[i]*zeta[i];
    }
    SendBuffer(ainv[iu], gamma[iu], beta[iu], zeta[iu], 
      sum_beta_gamma, sum_beta_zeta, k, j, tblock);
  }

  SaveCoefficients(ainv, gamma, zeta, k, j, il, iu);
}

// delta -> x
// zeta -> r
// alpha -> c
// diagU -> a
// gamma -> v

template<typename T1, typename T2>
void ImplicitSolver::PeriodicBackwardSubstitution(
  std::vector<T1> &diagU,
  std::vector<T2> &delta, 
  int kl, int ku, int jl, int ju, int il, int iu)
{
  T2 delta_last;
  std::vector<T1> gamma(diagU.size()), ainv(diagU.size());
  std::vector<T2> zeta(diagU.size());

  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j) {
      LoadCoefficients(ainv, gamma, zeta, k, j, il, iu);
      if (last_block) {  // I'm the last block
        // a[n-1] = 0
        diagU[iu-1].setZero();
        //           -1
        // x[n] = c[n]*r[n]
        delta[iu]  = ainv[iu]*zeta[iu];
        delta_last = delta[iu];
      } else {
        RecvBuffer(delta[iu+1], delta_last, k, j, tblock);
        delta[iu] = ainv[iu]*(zeta[iu] - diagU[iu]*delta[iu+1] - gamma[iu]*delta_last);
      }

      // backward substitution
      for (int i = iu-1; i >= il; --i)
        //           -1
        // x[i] = c[i]*(r[i] - a[i]*x[i+1] - v[i]*x[n])
        delta[i] = ainv[i]*(zeta[i] - diagU[i]*delta[i+1] - gamma[i]*delta_last);

      // update conserved variables
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

      if (!first_block)
        SendBuffer(delta[il], delta_last, k, j, bblock);
    }

#ifdef MPI_PARALLEL
  MPI_Status status;

  if (!last_block)
    if (tblock.snb.rank != Globals::my_rank)
      for (int k = kl; k <= ku; ++k)
        for (int j = jl; j <= ju; ++j)
          MPI_Wait(&req_send_data6_[k][j], &status);
  else
    if (tblock.snb.rank != Globals::my_rank)
      MPI_Wait(&req_send_data1_[0][0], &status);

  if (!first_block)
    if (bblock.snb.rank != Globals::my_rank)
      for (int k = kl; k <= ku; ++k)
        for (int j = jl; j <= ju; ++j)
          MPI_Wait(&req_send_data2_[k][j], &status);
#endif
}

#endif
