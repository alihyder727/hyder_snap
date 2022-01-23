#ifndef INVERSION_HPP
#define INVERSION_HPP

// C/C++ headers
#include <string>

// Athena++ header
#include "../athena.hpp"
#include "mcmc.hpp"

class MeshBlock;
class ParameterInput;
class RadioObservation;

class Inversion {
public:
  // data
  std::string task;
  MeshBlock *pmy_block;
  RadioObservation *pradio;

  // functions
  Inversion(MeshBlock *pmb, ParameterInput *pin);
  ~Inversion();

  // MCMC functions
  void MakeMCMCOutputs(std::string fname);
  void MCMCStep();
  void ResetChain();

private:
  // mcmc variables
  mcmc_opts opts_;
  mcmc_recs recs_;
  
  bool mcmc_initialized_;
};

#endif
