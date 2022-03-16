#ifndef INVERSION_HPP
#define INVERSION_HPP

// C/C++ headers
#include <string>

// Athena++ header
#include "../athena.hpp"
#include "mcmc.hpp"

// Eigen header files
#include "../math/eigen335/Eigen/Core"

class MeshBlock;
class ParameterInput;

class Inversion {
public:
  // general data
	std::string task;
	Inversion *next, *prev;
  MeshBlock *pmy_block;

  Eigen::VectorXd target;
  Eigen::MatrixXd icov;

  int nwalker, ndim, nvalue;

  // derived class data
	std::vector<int> ix;
	std::vector<Real> plevel;

  // functions
  Inversion(MeshBlock *pmb, ParameterInput *pin);
  virtual ~Inversion();
  virtual void ReadObservationFile(std::string fname);
	virtual Real LogPosteriorProbability(Real const *par, Real *val, 
		int ndim_, int nvalue_, int kwalker) const = 0;
	virtual void UpdateAtmosphere(void *params, int k = -1, int j = -1) const = 0;

  // MCMC functions
	void InitializeChain(int nwalker, int ndim, int nvalue);
  void MakeMCMCOutputs(std::string fname);
  void MCMCStep();
  void ResetChain();

protected:
	// mcmc initial positions
  Real **init_pos_;

private:
  // mcmc variables
  mcmc_opts opts_;
  mcmc_recs recs_;
  bool mcmc_initialized_;
};

void calculate_fit_target(MeshBlock *pmb, Real *val, int nvalue,
    int k, int j, bool differential = false);

#endif
