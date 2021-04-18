#ifndef COMMUNICATION_HPP
#define COMMUNICATION_HPP

// MPI headers
#ifdef MPI_PARALLEL
  #include <mpi.h>
#endif

// C/C++ headers
#include <cstring>

// Athena++ headers
#include "../../globals.hpp"
#include "../../mesh/mesh.hpp"
#include "../hydro.hpp"
#include "implicit_solver.hpp"

template<typename T>
void ImplicitSolver::SendBuffer(T const& a, int k, int j, NeighborBlock nb) {
  int s1 = a.size();
  int phy = k << 10 | j << 3 | 0;

  memcpy(buffer_, a.data(), s1*sizeof(Real));

  if (nb.snb.rank != Globals::my_rank) { // MPI boundary
#ifdef MPI_PARALLEL
    int tag = CreateMPITag(nb.snb.gid, pmy_hydro->pmy_block->gid, phy);
    MPI_Isend(buffer_, s1, MPI_ATHENA_REAL, nb.snb.rank, tag, MPI_COMM_WORLD,
      &req_send_data1_[k][j]);
#endif
  } else { // local boundary
    MeshBlock *pbl = pmy_hydro->pmy_block->pmy_mesh->FindMeshBlock(bblock.snb.gid);
    if (mydir == X1DIR)
      std::memcpy(pbl->phydro->pimp1->buffer_, buffer_, s1*sizeof(Real));
    else if (mydir == X2DIR)
      std::memcpy(pbl->phydro->pimp2->buffer_, buffer_, s1*sizeof(Real));
    else
      std::memcpy(pbl->phydro->pimp3->buffer_, buffer_, s1*sizeof(Real));
  }
}

template<typename T1, typename T2>
void ImplicitSolver::SendBuffer(T1 const &a, T2 const &b, int k, int j, NeighborBlock nb) {
  int s1 = a.size(), s2 = b.size();
  int phy = k << 10 | j << 3 | 1;

  memcpy(buffer_, a.data(), s1*sizeof(Real));
  memcpy(buffer_ + s1, b.data(), s2*sizeof(Real));

  if (nb.snb.rank != Globals::my_rank) { // MPI boundary
#ifdef MPI_PARALLEL
    int tag = CreateMPITag(nb.snb.gid, pmy_hydro->pmy_block->gid, phy);
    MPI_Isend(buffer_, s1+s2, MPI_ATHENA_REAL, nb.snb.rank, tag, MPI_COMM_WORLD,
      &req_send_data2_[k][j]);
#endif
  } else { // local boundary
    MeshBlock *pbl = pmy_hydro->pmy_block->pmy_mesh->FindMeshBlock(bblock.snb.gid);
    if (mydir == X1DIR)
      std::memcpy(pbl->phydro->pimp1->buffer_, buffer_, (s1+s2)*sizeof(Real));
    else if (mydir == X2DIR)
      std::memcpy(pbl->phydro->pimp2->buffer_, buffer_, (s1+s2)*sizeof(Real));
    else
      std::memcpy(pbl->phydro->pimp3->buffer_, buffer_, (s1+s2)*sizeof(Real));
  }
}

template<typename T1, typename T2, typename T3, 
         typename T4, typename T5, typename T6, typename T7>
void ImplicitSolver::SendBuffer(T1 const& a, T2 const&b, T3 const& c, T4 const& d,
  T5 const& e, T6 const& f, T7 const& g, int k, int j, NeighborBlock nb) {
  int s1 = a.size(), s2 = b.size(), s3 = c.size(), s4 = d.size();
  int s5 = e.size(), s6 = f.size(), s7 = g.size();
  int phy = k << 10 | j << 3 | 5;

  Real *it = buffer_;
  memcpy(it, a.data(), s1*sizeof(Real));
  it += s1;
  memcpy(it, b.data(), s2*sizeof(Real));
  it += s2;
  memcpy(it, c.data(), s3*sizeof(Real));
  it += s3;
  memcpy(it, d.data(), s4*sizeof(Real));
  it += s4;
  memcpy(it, e.data(), s5*sizeof(Real));
  it += s5;
  memcpy(it, f.data(), s6*sizeof(Real));
  it += s6;
  memcpy(it, g.data(), s7*sizeof(Real));

  int st = s1+s2+s3+s4+s5+s6+s7;

  if (nb.snb.rank != Globals::my_rank) { // MPI boundary
#ifdef MPI_PARALLEL
    int tag = CreateMPITag(nb.snb.gid, pmy_hydro->pmy_block->gid, phy);
    MPI_Isend(buffer_, st, MPI_ATHENA_REAL, nb.snb.rank, tag, MPI_COMM_WORLD,
      &req_send_data7_[k][j]);
#endif
  } else { // local boundary
    MeshBlock *pbl = pmy_hydro->pmy_block->pmy_mesh->FindMeshBlock(bblock.snb.gid);
    if (mydir == X1DIR)
      std::memcpy(pbl->phydro->pimp1->buffer_, buffer_, st*sizeof(Real));
    else if (mydir == X2DIR)
      std::memcpy(pbl->phydro->pimp2->buffer_, buffer_, st*sizeof(Real));
    else
      std::memcpy(pbl->phydro->pimp3->buffer_, buffer_, st*sizeof(Real));
  }
}

template<typename T>
void ImplicitSolver::RecvBuffer(T &a, int k, int j, NeighborBlock nb) {
  int s1 = a.size();
  int phy = k << 10 | j << 3 | 0;
#ifdef MPI_PARALLEL
  MPI_Status status;
#endif

  if (nb.snb.rank != Globals::my_rank) { // MPI boundary
#ifdef MPI_PARALLEL
    int tag = CreateMPITag(pmy_hydro->pmy_block->gid, nb.snb.gid, phy);
    MPI_Recv(buffer_, s1, MPI_ATHENA_REAL, nb.snb.rank, tag, MPI_COMM_WORLD, &status);
#endif
  } // local boundary

  memcpy(a.data(), buffer_, s1*sizeof(Real));
}


template<typename T1, typename T2>
void ImplicitSolver::RecvBuffer(T1 &a, T2 &b, int k, int j, NeighborBlock nb) {
  int s1 = a.size(), s2 = b.size();
  int phy = k << 10 | j << 3 | 1;
#ifdef MPI_PARALLEL
  MPI_Status status;
#endif

  if (nb.snb.rank != Globals::my_rank) { // MPI boundary
#ifdef MPI_PARALLEL
    int tag = CreateMPITag(pmy_hydro->pmy_block->gid, nb.snb.gid, phy);
    MPI_Recv(buffer_, s1+s2, MPI_ATHENA_REAL, nb.snb.rank, tag, MPI_COMM_WORLD, &status);
#endif
  } // local boundary

  memcpy(a.data(), buffer_, s1*sizeof(Real));
  memcpy(b.data(), buffer_ + s1, s2*sizeof(Real));
}

template<typename T1, typename T2, typename T3,
         typename T4, typename T5, typename T6, typename T7>
void ImplicitSolver::RecvBuffer(T1 &a, T2 &b, T3 &c, T4 &d, 
  T5 &e, T6 &f, T7 &g, int k, int j, NeighborBlock nb) {
  int s1 = a.size(), s2 = b.size(), s3 = c.size(), s4 = d.size();
  int s5 = e.size(), s6 = f.size(), s7 = g.size();
  int phy = k << 10 | j << 3 | 5;
#ifdef MPI_PARALLEL
  MPI_Status status;
#endif

  int st = s1+s2+s3+s4+s5+s6+s7;

  if (nb.snb.rank != Globals::my_rank) { // MPI boundary
#ifdef MPI_PARALLEL
    int tag = CreateMPITag(pmy_hydro->pmy_block->gid, nb.snb.gid, phy);
    MPI_Recv(buffer_, st, MPI_ATHENA_REAL, nb.snb.rank, tag, MPI_COMM_WORLD, &status);
#endif
  } // local boundary

  Real *it = buffer_;
  memcpy(a.data(), it, s1*sizeof(Real));
  it += s1;
  memcpy(b.data(), it, s2*sizeof(Real));
  it += s2;
  memcpy(c.data(), it, s3*sizeof(Real));
  it += s3;
  memcpy(d.data(), it, s4*sizeof(Real));
  it += s4;
  memcpy(e.data(), it, s5*sizeof(Real));
  it += s5;
  memcpy(f.data(), it, s6*sizeof(Real));
  it += s6;
  memcpy(g.data(), it, s7*sizeof(Real));
}

template<typename T1, typename T2>
void ImplicitSolver::SaveCoefficients(std::vector<T1> &a, std::vector<T2> &b,
  int k, int j, int il, int iu) {
  for (int i = il; i <= iu; ++i) {
    int s1 = a[i].size(), s2 = b[i].size();
    memcpy(coefficients_[k][j][i], a[i].data(), s1*sizeof(Real));
    memcpy(coefficients_[k][j][i] + s1, b[i].data(), s2*sizeof(Real));
  }
}

template<typename T1, typename T2, typename T3>
void ImplicitSolver::SaveCoefficients(std::vector<T1> &a, std::vector<T2> &b,
  std::vector<T3> &c, int k, int j, int il, int iu) {
  for (int i = il; i <= iu; ++i) {
    int s1 = a[i].size(), s2 = b[i].size(), s3 = c[i].size();
    memcpy(coefficients_[k][j][i], a[i].data(), s1*sizeof(Real));
    memcpy(coefficients_[k][j][i]+s1, b[i].data(), s2*sizeof(Real));
    memcpy(coefficients_[k][j][i]+s1+s2, c[i].data(), s3*sizeof(Real));
  }
}

template<typename T1, typename T2>
void ImplicitSolver::LoadCoefficients(std::vector<T1> &a, std::vector<T2> &b,
  int k, int j, int il, int iu) {
  for (int i = il; i <= iu; ++i) {
    int s1 = a[i].size(), s2 = b[i].size();
    memcpy(a[i].data(), coefficients_[k][j][i], s1*sizeof(Real));
    memcpy(b[i].data(), coefficients_[k][j][i] + s1, s2*sizeof(Real));
  }
}

template<typename T1, typename T2, typename T3>
void ImplicitSolver::LoadCoefficients(std::vector<T1> &a, std::vector<T2> &b,
  std::vector<T3> &c, int k, int j, int il, int iu) {
  for (int i = il; i <= iu; ++i) {
    int s1 = a[i].size(), s2 = b[i].size(), s3 = c[i].size();
    memcpy(a[i].data(), coefficients_[k][j][i], s1*sizeof(Real));
    memcpy(b[i].data(), coefficients_[k][j][i]+s1, s2*sizeof(Real));
    memcpy(c[i].data(), coefficients_[k][j][i]+s1+s2, s3*sizeof(Real));
  }
}

#endif
