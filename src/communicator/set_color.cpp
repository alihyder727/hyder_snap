/** @file set_color.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Thursday Apr 14, 2022 11:41:57 EDT
 * @bug No known bugs.
 */

#include "communicator.hpp"
#include "../globals.hpp"
#include "../mesh/mesh.hpp"

void Communicator::setColor(CoordinateDirection dir) {
  MeshBlock *pmb = pmy_block_;
  NeighborBlock bblock, tblock;
  pmb->FindNeighbors(dir, bblock, tblock);

  if (dir == X1DIR) {
    if (pmb->block_size.x1min <= pmb->pmy_mesh->mesh_size.x1min) {
      bblock.snb.gid = -1;
      bblock.snb.rank = -1;
    } 
    if (pmb->block_size.x1max >= pmb->pmy_mesh->mesh_size.x1max) {
      tblock.snb.gid = -1;
      tblock.snb.rank = -1;
    }
  } else if (dir == X2DIR) {
    if (pmb->block_size.x2min <= pmb->pmy_mesh->mesh_size.x2min) {
      bblock.snb.gid = -1;
      bblock.snb.rank = -1;
    }
    if (pmb->block_size.x2max >= pmb->pmy_mesh->mesh_size.x2max) {
      tblock.snb.gid = -1;
      tblock.snb.rank = -1;
    }
  } else { // X3DIR
    if (pmb->block_size.x3min <= pmb->pmy_mesh->mesh_size.x3min) {
      bblock.snb.gid = -1;
      bblock.snb.rank = -1;
    }
    if (pmb->block_size.x3max >= pmb->pmy_mesh->mesh_size.x3max) {
      tblock.snb.gid = -1;
      tblock.snb.rank = -1;
    }
  }

#ifdef MPI_PARALLEL
  MPI_Allgather(&bblock.snb.rank, 1, MPI_INT, brank_, 1, MPI_INT, MPI_COMM_WORLD);
#else
  brank_[0] = -1;
#endif

  std::fill(color_, color_ + Globals::nranks, -1);
  int c = 0;
  for (int i = 0; i < Globals::nranks; ++i) {
    //color[i] = brank_[i] == -1 ? color[i] : color[brank_[i]];
    if (brank_[i] == -1) {
      if (color_[i] == -1)
        color_[i] = c++;
    } else
      color_[i] = color_[brank_[i]];
  }

#ifdef MPI_PARALLEL
  if (comm_ != MPI_COMM_WORLD)
    MPI_Comm_free(&comm_);
  MPI_Comm_split(MPI_COMM_WORLD, color_[Globals::my_rank], Globals::my_rank, &comm_);
#endif
}
