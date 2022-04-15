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
  pmb->pdebug->Leave();
}

Communicator::~Communicator()
{
  delete[] color_;
  delete[] brank_;
}

int Communicator::getRank(CoordinateDirection dir)
{
  int r = 0;
  int b = brank_[Globals::my_rank];
  while (b != -1) {
    r++;
    b = brank_[b];
  }
  return r;
}

void Communicator::gatherData(Real *send, Real *recv, int size, CoordinateDirection dir)
{
#ifdef MPI_PARALLEL
  MPI_Comm comm;
  std::fill(color_, color_ + Globals::nranks, -1);
  MPI_Comm_split(MPI_COMM_WORLD, color_[Globals::my_rank], Globals::my_rank, &comm);
  MPI_Allgather(send, size, MPI_ATHENA_REAL, recv, size, MPI_ATHENA_REAL, comm);
  MPI_Comm_free(&comm);
#else
  memcpy(recv, send, size*sizeof(Real));
#endif
}
