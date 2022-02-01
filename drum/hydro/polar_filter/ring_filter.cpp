/** @file ring_filter.cpp
 * @brief implements a ring filter for polar region
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Sunday Aug 01, 2021 12:58:07 EDT
 * @bug No known bugs.
 */

// C++ headers
#include <iostream>
#include <sstream>

// Athena++ headers
#include "../../mesh/mesh.hpp"
#include "../../globals.hpp"
#include "../../debugger/debugger.hpp"
#include "ring_filter.hpp"

RingFilter::RingFilter(Hydro *phydro):
  pmy_hydro(phydro), my_rank(0)
{
  phydro->pmy_block->pdebug->Enter("RingFilter");
  std::stringstream msg;
  MeshBlock *pmb = phydro->pmy_block;
  Mesh *pm = pmb->pmy_mesh;
  int nc1 = pmb->block_size.nx1, nc3 = pmb->block_size.nx3;
  int nx3 = pm->mesh_size.nx3;

  if (strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    if (ceil(log2(nx3)) != floor(log2(nx3))) {
      msg << "### FATAL ERROR in RingFilter::RingFilter" << std::endl
          << "Number of longitudinal cells must be a power of 2";
      ATHENA_ERROR(msg);
    }
  }

  nlevel = int(log2(nx3/8));
  lrank_ = new int [Globals::nranks];
  color_ = new int [Globals::nranks];
  buffer_send_ = new Real [NHYDRO*nc1*nlevel*nc3];
  buffer_recv_ = new Real [NHYDRO*nc1*nlevel*nx3];

  hydro_mean.NewAthenaArray(NHYDRO,nx3,nlevel,nc1);
  phydro->pmy_block->pdebug->Leave();
}

RingFilter::~RingFilter()
{
  delete[] lrank_;
  delete[] color_;
  delete[] buffer_send_;
  delete[] buffer_recv_;
}

void RingFilter::FindNeighbors()
{
  NeighborBlock lblock, rblock;
  lblock.snb.gid = -1;
  lblock.snb.rank = -1;

  rblock.snb.gid = -1;
  rblock.snb.rank = -1;

  MeshBlock *pmb = pmy_hydro->pmy_block;
  for (int n = 0; n < pmb->pbval->nneighbor; ++n) {
    NeighborBlock& nb = pmb->pbval->neighbor[n];
    if ((nb.ni.ox1 == 0) && (nb.ni.ox2 == 0) && (nb.ni.ox3 == -1)) {
      lblock = nb;
    } if ((nb.ni.ox1 == 0) && (nb.ni.ox2 == 0) && (nb.ni.ox3 == 1)) {
      rblock = nb;
    }
  }

  // set first block
  if (pmb->block_size.x3min == pmb->pmy_mesh->mesh_size.x3min) {
    lblock.snb.gid = -1;
    lblock.snb.rank = -1;
  }

#ifdef MPI_PARALLEL
  MPI_Allgather(&lblock.snb.rank, 1, MPI_INT, lrank_, 1, MPI_INT, MPI_COMM_WORLD);
#else
  lrank_[0] = -1;
#endif

  int c = 0;
  for (int i = 0; i < Globals::nranks; ++i) {
    if (lrank_[i] == -1)
      color_[i] = c++;
    else
      color_[i] = color_[lrank_[i]];
  }
}

void RingFilter::PopulateConserved(AthenaArray<Real> const& u, bool north_pole) {
  // std::cout << "I'm in PopulateConserved" << std::endl;
  // pack data into 1D array
  MeshBlock *pmb = pmy_hydro->pmy_block;
  Mesh *pm = pmb->pmy_mesh;
  int nc1 = pmb->block_size.nx1, nc3 = pmb->block_size.nx3;
  int nx3 = pm->mesh_size.nx3;
  int ssize = NHYDRO*nc1*nlevel*nc3;
  int ks = pmb->ks, js = pmb->js, is = pmb->is;
  int ke = pmb->ke, je = pmb->je, ie = pmb->ie;
  int p = 0;

  if (north_pole) {
    for (int n = 0; n < NHYDRO; ++n)
      for (int k = ks; k <= ke; ++k)
        for (int j = js; j < js + nlevel; ++j) 
          for (int i = is; i <= ie; ++i)
            buffer_send_[p++] = u(n,k,j,i);
  } else {  // south pole
    for (int n = 0; n < NHYDRO; ++n)
      for (int k = ks; k <= ke; ++k)
        for (int j = je; j > je - nlevel; --j)
          for (int i = is; i <= ie; ++i)
            buffer_send_[p++] = u(n,k,j,i);
  }

  // gather all data on the latitude ring
#ifdef MPI_PARALLEL
  MPI_Comm comm_ring;
  MPI_Comm_split(MPI_COMM_WORLD, color_[Globals::my_rank], Globals::my_rank, &comm_ring);
  // assuming correct ordering
  MPI_Allgather(buffer_send_, ssize, MPI_ATHENA_REAL, buffer_recv_, ssize, 
        MPI_ATHENA_REAL, comm_ring);
  MPI_Comm_rank(comm_ring, &my_rank);
  MPI_Comm_free(&comm_ring);
#endif

  // unpack data
  p = 0;
  for (int r = 0; r < nx3/nc3; ++r)
    for (int n = 0; n < NHYDRO; ++n)
      for (int k = 0; k < nc3; ++k)
        for (int j = 0; j < nlevel; ++j)
          for (int i = 0; i < nc1; ++i)
            hydro_mean(n,r*nc3 + k,j,i) = buffer_recv_[p++];

  /*if (Globals::my_rank == 0)
    for (int k = 0; k < nx3; ++k)
      std::cout << hydro_mean(IV1,k,0,0) << std::endl;*/
}

void RingFilter::ApplyPolarFilter(AthenaArray<Real> &u)
{
  if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") != 0)
    return;

  MeshBlock *pmb = pmy_hydro->pmy_block;
  FindNeighbors();
  PopulateConserved(u, 1);
  // south pole
  if (pmb->block_size.x2min < 0.01)
    ApplyRingFilter(u, 1);

  // north pole
  PopulateConserved(u, 0);
  if (pmb->block_size.x2max > 3.13)
    ApplyRingFilter(u, 0);
}
