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

RadioData::RadioData(MeshBlock *pmb, ParameterInput *pin,
    Real grav_, int iH2O_, int iNH3_,
    Real *z1_, Real *p1_, Real *t1_, int nx1_, bool test_):
  grav(grav_), iH2O(iH2O_), iNH3(iNH3_), nx1(nx1_), testrun_(test_)
{
  std::stringstream msg;
  is = pmb->is; js = pmb->js; ks = pmb->ks;
  ie = pmb->ie; je = pmb->je; ke = pmb->ke;

  pthermo = pmb->pthermo;
  phydro = pmb->phydro;
  prad = pmb->prad;
  pcoord = pmb->pcoord;

  std::vector<Direction> out_dir = prad->GetOutgoingRays();
  int nwave = prad->GetNumBands();
  int nangle = out_dir.size();

  NewCArray(target, nwave, 3);
  NewCArray(icov, nwave, 3, 3);
  std::fill(**icov, **icov + nwave*9, 0.);

  std::string obsfile = pin->GetOrAddString("problem", "obsfile", "");
  Real delta = pin->GetReal("problem", "delta");
  Real eps = pin->GetReal("problem", "eps");

  // solve coefficients using least square
  Eigen::MatrixXd B(nangle,3);
  Eigen::VectorXd b(nangle), sol(3);
  for (int i = 0; i < nangle; ++i) {
    B(i,0) = 1.;
    B(i,1) = 1. - out_dir[i].mu;
    B(i,2) = sqr(1. - out_dir[i].mu);
  }

  if (obsfile != "") {
    ReadObservationFile(obsfile.c_str());
  } else {
    RadiationBand *pband = prad->pband;
    int i = 0;
    while (pband != NULL) {
      for (int j = 0; j < nangle; ++j)
        b(j) = pband->btoa(j,ks,js);
      sol = B.colPivHouseholderQr().solve(b);

      icov[i][0][0] = sqr((sol(0)-300.)*delta) + sqr(eps);
      icov[i][1][1] = sqr(eps);
      icov[i][2][2] = sqr(2*eps);

      if (i == 0) for (int j = 0; j < nwave*9; ++j) (**icov)[j] *= 4;

      Eigen::Map<Eigen::MatrixXd> cov(*icov[i], 3, 3);
      Eigen::LLT<Eigen::MatrixXd> LLTofCov(cov);
      if (LLTofCov.info() == Eigen::NumericalIssue)
        throw std::runtime_error("Covariance matrix is not semi-positive definitie!");
      cov = cov.inverse();

      for (int j = 0; j < 3; ++j)
        target[i][j] = sol(j);
      i++;
      pband = pband->next;
    }
  }

  pdiv = Vectorize<Real>(pin->GetString("problem", "pdiv").c_str());
  zfrac = Vectorize<Real>(pin->GetString("problem", "zfrac").c_str());
  int nsample = zfrac.size();

  zdiv.resize(nsample+1);
  for (int i = 0; i < nsample; ++i)
    zdiv[i] = interp1(pdiv[i]*1.E5, z1_, p1_, nx1_);  // bar to pa
  zdiv[nsample] = pcoord->x1v(ie);

  z1 = new Real [nx1];
  p1 = new Real [nx1];
  t1 = new Real [nx1];
  memcpy(z1, z1_, nx1_*sizeof(Real));
  memcpy(p1, p1_, nx1_*sizeof(Real));
  memcpy(t1, t1_, nx1_*sizeof(Real));

  // temperature correlation 
  Tstd = pin->GetReal("problem", "Tstd");
  Tlen = pin->GetReal("problem", "Tlen")*1.E3; // km -> m

  // temperature perturbation
  TpSample = Vectorize<Real>(pin->GetString("problem", "Tp").c_str());

  // check size
  if (TpSample.size() != nsample) {
    msg << "### FATAL ERROR in RadioData" << std::endl 
        << "size of temperature perturbation (Tp) should be "
        << nsample;
    throw std::runtime_error(msg.str().c_str());
  }

  // ammonia correlation
  NH3std = pin->GetReal("problem", "NH3std")/1.E3;  // g/kg -> kg/kg
  NH3len = pin->GetReal("problem", "NH3len")*1.E3; // km -> m

  // read ammonia perturbation
  NH3pSample = Vectorize<Real>(pin->GetString("problem", "NH3p").c_str());

  for (int i = 0; i < NH3pSample.size(); ++i)
    NH3pSample[i] /= 1.E3;  // g/kg -> kg/kg

  // check size
  if (NH3pSample.size() != nsample) {
    msg << "### FATAL ERROR in RadioData" << std::endl
        << "size of ammonia perturbation (NH3p) should be "
        << nsample;
    throw std::runtime_error(msg.str().c_str());
  }

}

RadioData::~RadioData()
{
  FreeCArray(target);
  FreeCArray(icov);
  delete[] z1;
  delete[] p1;
  delete[] t1;
}

#define MAX_LINE 512
void RadioData::WriteObservationFile(char const *fname)
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
}

void RadioData::ReadObservationFile(char const *fname)
{
  std::stringstream msg;
  FILE *fp = fopen(fname, "r");
  if (fp == NULL) {
    msg << "### FATAL ERROR in RadioData::ReadObseravtionFile" << std::endl 
        << fname << " cannot be opened.";
    throw std::runtime_error(msg.str().c_str());
  }
  char line[MAX_LINE], *pl;

  int chunks, rows, cols;
  pl = NextLine(line, MAX_LINE, fp);

  int nwave = prad->GetNumBands();

  // target values
  sscanf(pl, "%d%d", &rows, &cols);
  assert(rows == nwave && cols == 3);
  for (int i = 0; i < rows; ++i) {
    pl = NextLine(line, MAX_LINE, fp);
    char *p = strtok(pl, " ");
    for (int j = 0; j < cols; ++j) {
      sscanf(p, "%lf", &target[i][j]);
      p = strtok(NULL, " ");
    }
  }

  // covariance matrix
  pl = NextLine(line, MAX_LINE, fp);
  sscanf(pl, "%d%d%d", &chunks, &rows, &cols);
  assert(chunks == nwave && rows == 3 && cols == 3);
  for (int i = 0; i < chunks; ++i)
    for (int j = 0; j < rows; ++j) {
      pl = NextLine(line, MAX_LINE, fp);
      char *p = strtok(pl, " ");
      for (int k = 0; k < cols; ++k) {
        sscanf(p, "%lf", &icov[i][j][k]);
        p = strtok(NULL, " ");
      }
    }

  // additional TP data to fit
  pl = NextLine(line, MAX_LINE, fp);
  if (strlen(pl) > 0) {
    sscanf(pl, "%d%d", &rows, &cols);
    tpdata.resize(rows);
    assert(cols == 3);
    for (int i = 0; i < rows; ++i) {
      pl = NextLine(line, MAX_LINE, fp);
      sscanf(pl, "%lf%lf%lf", &tpdata[i].P, &tpdata[i].T, &tpdata[i].ERR);
      tpdata[i].P *= 1.E5;  // bar to pa
    }
  }

  fclose(fp);
}
#undef MAX_LINE
