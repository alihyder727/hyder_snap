// Athena++ header files
#include "chemistry_solver.hpp"

template<typename T>
void ChemistryBase<T>::IntegrateDense(AthenaArray<Real> &u, AthenaArray<Real> &c,
  Real time, Real dt)
{
  MeshBlock *pmb = pmy_block;
  Thermodynamics *pthermo = pmb->pthermo;
  int ks = pmb->ks, js = pmb->js, is = pmb->is;
  int ke = pmb->ke, je = pmb->je, ie = pmb->ie;

  int size = u.GetDim4() + c.GetDim4();

  Real* q = new Real [size];
  Real* q1 = new Real [size];
  Real* q2 = new Real [size];
  Real mu_d = Thermodynamics::Rgas/pthermo->GetRd();

  T* pchem = static_cast<T*>(this);
  Eigen::Matrix<Real, T::Solver::Size, T::Solver::Size> Jac;
  Eigen::Matrix<Real, T::Solver::Size, 1> Rate, Sol;

  SetTotalCv(u, pmb->ppart, is, ie, js, je, ks, ke);

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Rate.setZero();
        Jac.setZero();

        pthermo->ConservedToChemical(q, u.at(k,j,i));
        mol_(k,j,i) += q[IPR]/(Thermodynamics::Rgas*q[IDN]);

        Real qd = 1.;
        for (int n = 1; n <= NVAPOR; ++n) qd -= q[n];
        for (int n = 0; n < c.GetDim4(); ++n) {
          Real epsc = pmy_part->GetMassRatio(n, mu_d);
          q[NHYDRO+n] = (c(n,k,j,i)/epsc)/u(0,k,j,i)*qd;
        }

        // make two copies
        std::memcpy(q1, q, size*sizeof(Real));
        std::memcpy(q2, q, size*sizeof(Real));

        pchem->AssembleReactionMatrix(Rate, Jac, q, cvt_(k,j,i)/mol_(k,j,i), time);
        // BDF1 solver
        Sol = pchem->solver.BDF1(q1, Rate, Jac, time, dt);
        for (int i = 0; i < T::Solver::Size; ++i)
          q1[qindex_[i]] += Sol(i);

        // TR-BDF2 solver
        Sol = pchem->solver.TRBDF2(q2, Rate, Jac, time, dt);
        for (int i = 0; i < T::Solver::Size; ++i)
          q2[qindex_[i]] += Sol(i);

        //pchem->ApplyChemicalLimits(q2, q0, cvt_(k,j,i)/mol_(k,j,i));

        // Blend solutions
        Real alpha = 1.;
        for (int i = 0; i < T::Solver::Size; ++i)
          if (q2[qindex_[i]] < 0.)
            alpha = std::min(alpha, q1[qindex_[i]]/(q1[qindex_[i]] - q2[qindex_[i]]));
        for (int i = 0; i < T::Solver::Size; ++i)
          q2[qindex_[i]] = (1. - alpha)*q1[qindex_[i]] + alpha*q2[qindex_[i]];

        //std::cout << "==== iter ends ===="<< std::endl;
        //std::cout << "iter = " << iter << std::endl;
        //std::cout << "norm = " << norm << std::endl;

        /* debug
        for (int n = 0; n < NHYDRO; ++n)
          std::cout << u(n,k,j,i) << " ";
        std::cout << c(0,k,j,i) << " ";
        std::cout << c(1,k,j,i) << std::endl;;
        std::cout << "i = " << i << " ===================" << std::endl;
        for (int n = 0; n < size; ++n)
          std::cout << q[n] << " ";
        std::cout << std::endl;
        for (int n = 0; n < size; ++n)
          std::cout << q1[n] << " ";
        std::cout << std::endl;
        std::cout << "=========================" << std::endl;*/

        // from molar mixing ratio to density
        pthermo->ChemicalToConserved(u.at(k,j,i), q2);
        for (int n = 0; n < c.GetDim4(); ++n) {
          Real epsc = pmy_part->GetMassRatio(n, mu_d);
          c(n,k,j,i) = (q2[NHYDRO+n]*epsc)/qd*u(0,k,j,i);
        }
        //for (int n = 0; n < NHYDRO; ++n)
        //  std::cout << u(n,k,j,i) << " ";
        //std::cout << c(0,k,j,i) << " ";
        //std::cout << c(1,k,j,i) << std::endl;;
        //std::cout << std::endl;
      }
  delete[] q;
  delete[] q1;
  delete[] q2;
}
