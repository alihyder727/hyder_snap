/** @file meshblock_additions.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @created Tuesday May 18, 2021 09:26:28 PDT
 * @bug No known bugs.
 */

#include "mesh.hpp"
#include "../coordinates/coordinates.hpp"

void MeshBlock::FindNeighbors(CoordinateDirection dir,
	NeighborBlock &bblock, NeighborBlock &tblock) {
  // set void bblock
	bblock.snb.gid = -1;
	bblock.snb.rank = -1;

	// set void tblock
	tblock.snb.gid = -1;
	tblock.snb.rank = -1;

  for (int n = 0; n < pbval->nneighbor; ++n) {
    NeighborBlock& nb = pbval->neighbor[n];
    if (dir == X1DIR) {
      if ((nb.ni.ox1 == -1) && (nb.ni.ox2 == 0) && (nb.ni.ox3 == 0))
        bblock = nb;
      if ((nb.ni.ox1 == 1) && (nb.ni.ox2 == 0) && (nb.ni.ox3 == 0))
        tblock = nb;
    } else if (dir == X2DIR) {
      if ((nb.ni.ox1 == 0) && (nb.ni.ox2 == -1) && (nb.ni.ox3 == 0))
        bblock = nb;
      if ((nb.ni.ox1 == 0) && (nb.ni.ox2 == 1) && (nb.ni.ox3 == 0))
        tblock = nb;
    } else { // X3DIR
      if ((nb.ni.ox1 == 0) && (nb.ni.ox2 == 0) && (nb.ni.ox3 == -1))
        bblock = nb;
      if ((nb.ni.ox1 == 0) && (nb.ni.ox2 == 0) && (nb.ni.ox3 == 1))
        tblock = nb;
    }
  }
}

