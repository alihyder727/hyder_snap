#ifndef INVERSION_HPP
#define INVERSION_HPP

// C/C++ headers
#include <string>

// Athena++ header
#include "../athena.hpp"
#include "mcmc.hpp"

class MeshBlock;
class ParameterInput;
template<typename T> class AthenaArray;

class Inversion {
public:
  // data
  std::string method;
  void *obj;  // objective values to fit

  // functions
  Inversion(MeshBlock *pmb, ParameterInput *pin);
  ~Inversion();
  void Initialize(Real **pos, int nwalker, int ndim, int nvalue);

  // MCMC functions
  void EnrollObjectives(ObjectiveFunction_t lnprob, void *myobj);
  void MakeMCMCOutputs(std::string fname);
  void MCMCStep(void *myobj);
  void ResetChain();

private:
  MeshBlock *pmy_block_;
  // objective function
  ObjectiveFunction_t lnprob_;

  // mcmc variables
  mcmc_opts opts_;
  mcmc_recs recs_;
  
  //bool objective_function_enrolled_;
  bool initialized_;
};

#endif
