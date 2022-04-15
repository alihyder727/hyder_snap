/** @file communicator.hpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Thursday Apr 14, 2022 11:39:35 EDT
 * @bug No known bugs.
 */

#ifndef COMMUNICATOR_HPP
#define COMMUNICATOR_HPP

// Athena++ header
#include "../athena.hpp"

class MeshBlock;

class Communicator {
public:
  Communicator(MeshBlock *pmb);
  ~Communicator();
  int getRank(CoordinateDirection dir);
  void setColor(CoordinateDirection dir);
  //void reduceData23(Real *send, Real *recv);
  void gatherData(Real *send, Real *recv, int size, CoordinateDirection dir);
  void gatherDataInPlace(Real *recv, int size, CoordinateDirection dir);

private:
  MeshBlock *pmy_block_;
  int *color_;
  int *brank_;
};

#endif /* end of include guard COMMUNICATOR_HPP */
