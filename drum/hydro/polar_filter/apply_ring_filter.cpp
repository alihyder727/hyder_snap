/** @file apply_ring_filter.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Sunday Aug 01, 2021 16:51:15 EDT
 * @bug No known bugs.
 */

// C/C++ header
#include <iostream>
#include <sstream>

// Athena++ headers
#include "ring_filter.hpp"
#include "../../mesh/mesh.hpp"
#include "../../globals.hpp"
#include "../../reconstruct/interpolation.hpp"

void RingFilter::ApplyRingFilter(AthenaArray<Real> &u, bool north_pole) {
  std::stringstream msg;
  MeshBlock *pmb = pmy_hydro->pmy_block;
  Mesh *pm = pmb->pmy_mesh;
  int nc1 = pmb->block_size.nx1, nc3 = pmb->block_size.nx3;
  int nx3 = pm->mesh_size.nx3;
  int ks = pmb->ks, js = pmb->js, is = pmb->is;
  int ke = pmb->ke, je = pmb->je, ie = pmb->ie;

  //std::cout << "I'm in ring filter" << std::endl;

  for (int j = 0; j < nlevel; ++j) {
    //if (Globals::my_rank == 0)
    //  std::cout << "I'm in level " << j << std::endl;
    int nchunk = 8 << j;
    int len = nx3/nchunk;
    //if (Globals::my_rank == 0) {
    //  std::cout << "nchunk = " << nchunk << std::endl;
    //  std::cout << "len = " << len << std::endl;
    //}

    for (int n = 0; n < NHYDRO; ++n)
      for (int i = 0; i < nc1; ++i) {
        // 1. averaging
        //std::cout << "hydro_mean = ";
        for (int r = 0; r < nchunk; ++r) {
          Real sum = 0.;
          for (int k = 0; k < len; ++k)
            sum += hydro_mean(n,r*len + k,j,i);
          hydro_mean(n,r*len,j,i) = sum/len;
          //std::cout << hydro_mean(n,r*len,j,i) << " ";
        }
        //std::cout << std::endl;

        // 2. reconstruction
        for (int r = 0; r < nchunk; ++r) {
          if (r*len > (1+my_rank)*nc3 || (r+1)*len < my_rank*nc3)
            continue;
          int r1 = (r-2+nchunk) % nchunk;
          int r2 = (r-1+nchunk) % nchunk;
          int r3 = (r) % nchunk;
          int r4 = (r+1) % nchunk;
          int r5 = (r+2) % nchunk;
          //std::cout << r1 << " " << r2 << " " << r3 << " " << r4 << " " << r5 << std::endl;
          Real ul = interp_weno5(hydro_mean(n,r1*len,j,i), hydro_mean(n,r2*len,j,i),
            hydro_mean(n,r3*len,j,i), hydro_mean(n,r4*len,j,i), hydro_mean(n,r5*len,j,i));
          Real ur = interp_weno5(hydro_mean(n,r5*len,j,i), hydro_mean(n,r4*len,j,i),
            hydro_mean(n,r3*len,j,i), hydro_mean(n,r2*len,j,i), hydro_mean(n,r1*len,j,i));
          Real um = hydro_mean(n,r3*len,j,i);
          //std::cout << ul << " " << ur << " " << um << std::endl;

          // 2.1 construct parabola
          Real A = 3.*(ul + ur - 2.*um);
          Real B = 2.*(3.*um - ur - 2.*ul);
          Real C = ul;

          // 2.2 interpolate
          if (north_pole) {
            for (int k = 1; k <= len; ++k) {
              int k_chunk = r*len + k - 1;
              int k_begin = my_rank*nc3;
              int k_end = (my_rank+1)*nc3;
              if (k_chunk >= k_begin && k_chunk < k_end)
                u(n,ks+k_chunk-k_begin,js+j,is+i) = 
                  A/(3.*nchunk*nchunk)*(3.*k*k - 3*k +1) + B/(2.*nchunk)*(2.*k - 1) + C;
              //if (Globals::my_rank == 0)
              //  std::cout << k << " " << k_chunk << " " << k_begin << std::endl;
              //std::cout << u(n,ks+k_chunk-k_begin,js+j,is+i) << std::endl;
            }
          } else {  // south pole
            for (int k = 1; k <= len; ++k) {
              int k_chunk = r*len + k - 1;
              int k_begin = my_rank*nc3;
              int k_end = (my_rank+1)*nc3;
              if (k_chunk >= k_begin && k_chunk < k_end)
                u(n,ks+k_chunk-k_begin,je-j,is+i) = 
                  A/(3.*nchunk*nchunk)*(3.*k*k - 3*k +1) + B/(2.*nchunk)*(2.*k - 1) + C;
            }
          }
        }
      }

    //if (Globals::my_rank == 0)
    //  for (int k = ks; k <= ke; ++k)
    //    std::cout << u(IV1,k,js,is) << std::endl;
  }
}

