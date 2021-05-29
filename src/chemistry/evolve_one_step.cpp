// Athena++ header files
#include "../thermodynamics/thermodynamics.hpp"
#include "../hydro/hydro.hpp"
#include "chemistry.hpp"
#include "chemistry_solver.hpp"

void Chemistry::EvolveOneStep(AthenaArray<Real> &u, Real time, Real dt)
{
  /* return if there's no phase
  if (NPHASE*NVAPOR == 0) return;

  MeshBlock *pmb = pmy_block_;
  Thermodynamics *pthermo = pmb->pthermo;
  int ks = pmb->ks, js = pmb->js, is = pmb->is;
  int ke = pmb->ke, je = pmb->je, ie = pmb->ie;

  Real q1[NHYDRO], q0[NHYDRO], dq[NMASS];
  ChemistrySolver solver;
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        // calculate density
        Real rho = 0., rho_hat = 0.;
        for (int n = 0; n < NMASS; ++n) {
          rho += u(n,k,j,i);
          rho_hat += u(n,k,j,i)/pthermo->GetMassRatio(n);
          q1[n] = u(n,k,j,i)/pthermo->GetMassRatio(n);
        }
        for (int n = 0; n < NMASS; ++n)
          q1[n] /= rho_hat;

        // calculate internal energy
        Real KE = 0.5*(u(IM1,k,j,i)*u(IM1,k,j,i) 
                     + u(IM2,k,j,i)*u(IM2,k,j,i)
                     + u(IM3,k,j,i)*u(IM3,k,j,i))/rho;
        Real uhat = (u(IEN,k,j,i) - KE)/(pthermo->GetRd()*rho_hat);
        pthermo->UpdateTPConservingU(q1, rho, uhat);

        int iter = 0;
        for (int n = 0; n < NHYDRO; ++n) q0[n] = q1[n];

        // debug
        //std::cout << "==== i = " << i;;
        //std::cout << "==== iter begins ===="<< std::endl;
        //for (int n = 0; n < NMASS; ++n)
        //  std::cout << q1[n] << " ";
        //std::cout << std::endl;

        // backward Euler
        std::fill(dq, dq + NMASS, 0.);
        Real norm;
        while (iter++ < max_iter_) {
          AssembleReactionMatrix(r0_, r1_, q1, time);

          for (int n = 0; n < NMASS; ++n)
            dq[n] = r0_[n] - dq[n]/dt;
          for (int n1 = 0; n1 < NMASS; ++n1)
            for (int n2 = 0; n2 < NMASS; ++n2)
              r1_[n1][n2] = 1./dt*identity_[n1][n2] - r1_[n1][n2];

          solver.SolveIndependent(r1_, dq);

          norm = 0;
          for (int n = 1; n < NMASS; ++n)
            norm += dq[n]*dq[n];
          if (std::sqrt(norm) < ftol_) break;

          for (int n = 1; n < NMASS; ++n) {
            q1[n] += dq[n]; 
            dq[n] = q1[n] - q0[n];
          }
          pthermo->UpdateTPConservingU(q1, rho, uhat);
          //std::cout << "norm = " << norm << std::endl;
        }

        //std::cout << "==== iter ends ===="<< std::endl;
        //std::cout << "iter = " << iter << std::endl;
        //std::cout << "norm = " << norm << std::endl;

        // debug
        //for (int n = 0; n < NMASS; ++n)
        //  std::cout << q1[n] << " ";
        //std::cout << std::endl;
        //std::cout << "=========================" << std::endl;

        // from molar mixing ratio to density
        Real sum = 1.;
        for (int n = 1; n < NMASS; ++n)
          sum += q1[n]*(pthermo->GetMassRatio(n) - 1.);
        for (int n = 1; n < NMASS; ++n)
          u(n,k,j,i) = rho*q1[n]*pthermo->GetMassRatio(n)/sum;
      }*/
}

/* TODO: TR-BDF2
// first step
rx = r0_;
Real a1 = 2./gamma_;
Real a2 = (2. - gamma_)/(1. - gamma_);
Real a3 = (1. - gamma_)/gamma_;
while (iter++ < max_iter_) {
  //std::cout << "==== iter = " << iter << " ===="<< std::endl;
  AssembleReactionMatrix(r0_, r1_, q, time);
  SolveChemistry(dq, a1/dt*identity_ - r1_, -a1*dq/dt + r0_ + rx, q);
  //std::cout << dq << std::endl;
  for (int n = 1; n < NMASS; ++n) {
    q[n] += dq(n);
    dq(n) = q[n] - q0[n];
  }
  pthermo->UpdateTPConservingU(q, rho, uhat);
}

// second step
iter = 0;
Vector dqdt = dq/dt;
dq.setZero();
for (int n = 0; n < NHYDRO; ++n)
  q0[n] = q[n];
while (iter++ < max_iter_) {
  AssembleReactionMatrix(r0_, r1_, q, time);
  SolveChemistry(dq, a2/dt*identity_ - r1_, -a2*dq/dt + r0_ + a3*dqdt, q);
  for (int n = 1; n < NMASS; ++n) {
    q[n] += dq(n);
    dq(n) = q[n] - q0[n];
  }
  pthermo->UpdateTPConservingU(q, rho, uhat);
}*/
