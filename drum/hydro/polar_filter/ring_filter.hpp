/** @file ring_filter.hpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Sunday Aug 01, 2021 13:15:00 EDT
 * @bug No known bugs.
 */

#ifndef RING_FILTER_HPP
#define RING_FILTER_HPP

// Athena++ headers
#include "../hydro.hpp"

class RingFilter {
public:
// data
  Hydro *pmy_hydro;
  AthenaArray<Real> hydro_mean;
  int nlevel;
  int my_rank;
  bool has_north_pole, has_south_pole;

// functions
  RingFilter(Hydro *phydro);
  ~RingFilter();
  void FindNeighbors();
  void PopulateConserved(AthenaArray<Real> const& u);
  void ApplyRingFilter(AthenaArray<Real> &u);

private:
  int *lrank_, *color_;   // left block rank and color
  Real *buffer_send_, *buffer_recv_;
};


#endif /* end of include guard RING_FILTER_HPP */
