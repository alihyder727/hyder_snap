// C/C++ headers
#include <iomanip>
#include <iostream>
#include <cstring>
#include <cassert>

// Athena++ headers
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "../hydro/hydro.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "../radiation/radiation.hpp"
#include "../utils/utils.hpp"
#include "../math/interpolation.h"
#include "../math/core.h"
#include "../math/eigen335/Eigen/Core"
#include "../math/eigen335/Eigen/Dense"
#include "radio_data.hpp"

RadioData::RadioData(MeshBlock *pmb, ParameterInput *pin):
  pmy_block(pmb)
{
  ATHENA_LOG("RadioData");
  std::string obsfile = pin->GetOrAddString("inversion", "obsfile", "");

  if (obsfile != "") {
    ReadObservationFile(obsfile.c_str());
    //target.setZero();
    //icov.setZero();
  } 
  std::cout << "- target: ";
  std::cout << target.transpose() << std::endl;
  std::cout << "- inverse covariance matrix" << std::endl;
  std::cout << icov << std::endl;

  // T correlation 
  Tstd = pin->GetReal("inversion", "Tstd");
  Tlen = pin->GetReal("inversion", "Tlen")*1.E3; // km -> m

  // X correlation
  Xstd = pin->GetReal("inversion", "Xstd")/1.E3;  // g/kg -> kg/kg
  Xlen = pin->GetReal("inversion", "Xlen")*1.E3; // km -> m

  // composition id
  ix = pin->GetInteger("inversion", "variables");
}

#define MAX_LINE 512
/*void RadioData::WriteObservationFile(char const *fname)
{
  FILE *fp = fopen(fname, "w");
  int nwave = prad->GetNumBands();

  // target values
  fprintf(fp, "# Target brightness temperatures:\n");
  fprintf(fp, "%10d%10d\n", nwave, 3);
  for (int i = 0; i < nwave; ++i) {
    for (int j = 0; j < 3; ++j)
      fprintf(fp, "%12.4f", target[i][j]);
    fprintf(fp, "\n");
  }

  // covariance matrix
  fprintf(fp, "# Inverse covariance matrix:\n");
  fprintf(fp, "%10d%10d%10d\n", nwave, 3, 3);
  for (int i = 0; i < nwave; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k)
        fprintf(fp, "%12.4f", icov[i][j][k]);
      fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
  }

  fprintf(fp, "# Additional TP data:\n");
  // additional TP data to fit
  if (tpdata.size() > 0) {
    fprintf(fp, "%10ld%10d\n", tpdata.size(), 3);
    for (int i = 0; i < tpdata.size(); ++i)
      fprintf(fp, "%12.4f%12.4f%12.4f\n", tpdata[i].P/1.E5, tpdata[i].T, tpdata[i].ERR);
  }

  fclose(fp);
}*/

void RadioData::ReadObservationFile(char const *fname)
{
  std::stringstream msg;
  FILE *fp = fopen(fname, "r");
  if (fp == NULL) {
    msg << "### FATAL ERROR in RadioData::ReadObseravtionFile" << std::endl 
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
