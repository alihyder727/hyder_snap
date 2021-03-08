#ifndef DEBUGGER_HPP
#define DEBUGGER_HPP

// C/C++ headers
#include <vector>
#include <string>

// Athena++ header
#include "../athena.hpp"

typedef int (*TestFunc_t)(Real);

class Debugger {
public:
// data
  MeshBlock *pmy_block;
  Debugger *prev, *next;

// functions
  Debugger(MeshBlock *pmb);
  ~Debugger();

  Debugger* StartTracking(std::string name);
  void Track3D(std::string name, TestFunc_t test, AthenaArray<Real>& var, int n);
  void Track1D(std::string name, TestFunc_t test, AthenaArray<Real>& var, int n, int k, int j);
  void DumpTracking(std::string name, int c1, int c2, int c3, char const* mode);

protected:
  std::string fname_;
  AthenaArray<Real> data_;
  std::vector<std::string> vnames_;
};

// test functions
int IsPositive(Real v);
int IsNumber(Real v);

#endif
