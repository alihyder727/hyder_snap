#ifndef JUNO_MWR_DATA_HPP_
#define JUNO_MWR_DATA_HPP_

// C/C++ headers
#include <vector>

// Athena++ headers
#include "../athena.hpp"

// Eigen header files
#include "../math/eigen335/Eigen/Core"

class MeshBlock;
class ParameterInput;

struct TPData {
  Real P, T, ERR;
};

class RadioData {
public:
  // data
  MeshBlock *pmy_block;

  Eigen::VectorXd target;
  Eigen::MatrixXd icov;

  Real Tstd, Tlen, Xstd, Xlen;
  int ix;
  //std::vector<TPData> tpdata;

  // functions
  RadioData(MeshBlock *pmb, ParameterInput *pin);
  void ReadObservationFile(char const *fname);
  //void WriteObservationFile(char const *fname);
};

void update_atm_profiles(MeshBlock *pmb,
    Real *PrSample, Real *TpSample, Real *XpSample, int nsample, int ix,
    Real Tstd, Real Tlen, Real Xstd, Real Xlen, Real P0, Real Z0 = 0.);

void calculate_fit_target(MeshBlock *pmb, int j, Real *val, int nvalue);

#endif
