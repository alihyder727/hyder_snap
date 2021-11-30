// C/C++ headers
#include <iomanip>
#include <iostream>
#include <cstring>
#include <cassert>

// Athena++ headers
#include "../mesh/mesh.hpp"
#include "../utils/utils.hpp"
//#include "../math/core.h"
#include "radio_observation.hpp"

RadioObservation::RadioObservation(Inversion *pinvt, ParameterInput *pin):
  pmy_invt_(pinvt)
{
  ATHENA_LOG("RadioObservation");
  std::string obsfile = pin->GetOrAddString("inversion", "obsfile", "none");

  if (obsfile != "none") {
    ReadObservationFile(obsfile.c_str());
    std::cout << "- target: ";
    std::cout << target.transpose() << std::endl;
    std::cout << "- inverse covariance matrix" << std::endl;
    std::cout << icov << std::endl;
  }

  // T correlation 
  Tstd_ = pin->GetReal("inversion", "Tstd");
  Tlen_ = pin->GetReal("inversion", "Tlen")*1.E3; // km -> m

  // X correlation
  Xstd_ = pin->GetReal("inversion", "Xstd")/1.E3;   // g/kg -> kg/kg
  Xlen_ = pin->GetReal("inversion", "Xlen")*1.E3;   // km -> m

  // power law coefficient
  chi_ = pin->GetOrAddReal("inversion", "chi", 0.1);

  // composition id
  ix = Vectorize<int>(pin->GetString("inversion", "Variables").c_str());

	// Pressure sample
	plevel = Vectorize<Real>(pin->GetString("inversion", "PrSample").c_str());
  std::cout << "- number of inversion variables: " << plevel.size()*ix.size() << std::endl;

	for (std::vector<Real>::iterator m = plevel.begin(); m != plevel.end(); ++m)
		(*m) *= 1.E5;	// bar -> pa

  // fit differential
  fit_differential_ = pin->GetOrAddBoolean("inversion", "differential", false);
}

#define MAX_LINE 512
void RadioObservation::ReadObservationFile(char const *fname)
{
  std::stringstream msg;
  FILE *fp = fopen(fname, "r");
  if (fp == NULL) {
    msg << "### FATAL ERROR in RadioObservation::ReadObseravtionFile" << std::endl 
        << fname << " cannot be opened.";
    ATHENA_ERROR(msg);
  }
  char line[MAX_LINE], *pl;

  int rows;
  // header
  pl = NextLine(line, MAX_LINE, fp);

  // target values
  sscanf(pl, "%d", &rows);
  target.resize(rows);
  icov.resize(rows, rows);

  for (int i = 0; i < rows; ++i) {
    pl = NextLine(line, MAX_LINE, fp);
    sscanf(pl, "%lf", &target(i));
  }

  // inverse covariance matrix
  for (int i = 0; i < rows; ++i) {
    pl = NextLine(line, MAX_LINE, fp);
    char *p = strtok(pl, " ");
    for (int j = 0; j < rows; ++j) {
      sscanf(p, "%lf", &icov(i,j));
      p = strtok(NULL, " ");
    }
  }

  fclose(fp);
}
#undef MAX_LINE
