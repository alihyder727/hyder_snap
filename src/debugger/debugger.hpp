#ifndef DEBUGGER_HPP
#define DEBUGGER_HPP

// C/C++ headers
#include <vector>
#include <string>
#include <sstream>

// Athena++ header
#include "../athena.hpp"
#include "../globals.hpp"

typedef int (*TestFunc_t)(Real);

class Debugger {
public:
// data
  MeshBlock *pmy_block;
  Debugger *prev, *next;
  std::stringstream msg;

// functions
  Debugger(MeshBlock *pmb);
  ~Debugger();

  Debugger* StartTracking(std::string name);
  void Track3D(std::string name, TestFunc_t test, AthenaArray<Real>& var, int n);
  void Track1D(std::string name, TestFunc_t test, AthenaArray<Real>& var, int n, int k, int j);
  void DumpTracking(std::string name, int c1, int c2, int c3, char const* mode);
  //void Enter(char const *name);
  void Enter(std::string name, std::string heil = "Initializing");
  void Call(std::string name) {
    Enter(name, "Calling");
  }
  void Leave();
  Debugger* WriteMessage(std::string str) const;
  

protected:
  std::string fname_;
  AthenaArray<Real> data_;
  std::vector<std::string> vnames_;
  std::vector<std::string> sections_;
  std::vector<std::string> idstack_next_;
};

void increment_id(std::string &str);

// test functions
int IsPositive(Real v);
int IsNumber(Real v);

#endif
