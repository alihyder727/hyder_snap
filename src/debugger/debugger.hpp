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

class MaterialPoint;

// DEBUG_LEVEL = 0 : no debug output
//             = 1 : output important steps
//             = 2 : output conservation check
//             = 3 : output detailed particle transfer and bounds check
class Debugger {
public:
// data
  static std::string const cgreen;
  static std::string const cend;

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
  //
  void Enter(std::string name, std::string heil = "Initializing");
  void Call(std::string name) {
    Enter(name, "Calling");
  }
  void Leave();
  void CheckConservation(std::string name, AthenaArray<Real> const& var,
      int is, int ie, int js, int je, int ks, int ke);
  void CheckParticleConservation(std::vector<std::string> const& cnames,
      std::vector<MaterialPoint> const& mp);
  //Debugger* WriteMessage(std::string str) const;
  

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
