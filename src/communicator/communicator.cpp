/** @file communicator.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Thursday Apr 14, 2022 11:48:11 EDT
 * @bug No known bugs.
 */

#include <iostream>

// Athena++ header
#include "communicator.hpp"
#include "../debugger/debugger.hpp"
#include "../mesh/mesh.hpp"
#include "../globals.hpp"

Communicator::Communicator(MeshBlock *pmb):
  pmy_block_(pmb)
{
  pmb->pdebug->Enter("Communicator");
  color_ = new int [Globals::nranks];
  brank_ = new int [Globals::nranks];
#ifdef MPI_PARALLEL
  comm_ = MPI_COMM_WORLD;
#endif
  pmb->pdebug->Leave();
}

Communicator::~Communicator()
{
  delete[] color_;
  delete[] brank_;

#ifdef MPI_PARALLEL
  if (comm_ != MPI_COMM_WORLD)
    MPI_Comm_free(&comm_);
#endif
}

int Communicator::getRank(CoordinateDirection dir) const
{
  int r = 0;
  int b = brank_[Globals::my_rank];
  while (b != -1) {
    r++;
    b = brank_[b];
  }
  return r;
}

void Communicator::gatherData(Real *send, Real *recv, int size) const
{
#ifdef MPI_PARALLEL
  MPI_Allgather(send, size, MPI_ATHENA_REAL, recv, size, MPI_ATHENA_REAL, comm_);
#else
  memcpy(recv, send, size*sizeof(Real));
#endif
}

//! \warning this function is unstable to some system.
void Communicator::gatherDataInPlace(Real *recv, int size) const
{
#ifdef MPI_PARALLEL
  MPI_Allgather(MPI_IN_PLACE, 0, 0, recv, size, MPI_ATHENA_REAL, comm_);
#endif
}

NeighborBlock const* Communicator::findBotNeighbor() const
{
  MeshBlock *pmb = pmy_block_;
  NeighborBlock *pbot = nullptr;
  for (int n = 0; n < pmb->pbval->nneighbor; ++n) {
    NeighborBlock* nb = pmb->pbval->neighbor + n;
    if ((nb->ni.ox1 == -1) && (nb->ni.ox2 == 0) && (nb->ni.ox3 == 0))
      pbot = nb;
  }
  return pbot;
}

NeighborBlock const* Communicator::findTopNeighbor() const
{
  MeshBlock *pmb = pmy_block_;
  NeighborBlock *ptop = nullptr;
  for (int n = 0; n < pmb->pbval->nneighbor; ++n) {
    NeighborBlock* nb = pmb->pbval->neighbor + n;
    if ((nb->ni.ox1 == 1) && (nb->ni.ox2 == 0) && (nb->ni.ox3 == 0))
      ptop = nb;
  }
  return ptop;
}
