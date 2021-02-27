//! \file pnetcdf.cpp
//  \brief writes output data in parallel-NETCDF format.
//  Data is written in RECTILINEAR_GRID geometry, in BINARY format, and in FLOAT type
//  Writes one file per timestep.

// C/C++ headers
#include <sstream>
#include <iostream>
#include <string>
#include <stdexcept>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>

// Athena++ headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../math/core.h"
#include "../coordinates/coordinates.hpp"
#include "outputs.hpp"

// Only proceed if PNETCDF output enabled
#ifdef PNETCDFOUTPUT

// External library headers
#include <mpi.h>
#include <pnetcdf.h>

#define ERR {if(err!=NC_NOERR){printf("Error at %s:%d : %s\n", __FILE__,__LINE__, ncmpi_strerror(err));}}

//----------------------------------------------------------------------------------------
// PnetcdfOutput constructor
// destructor - not needed for this derived class

PnetcdfOutput::PnetcdfOutput(OutputParameters oparams)
  : OutputType(oparams)
{
}

//----------------------------------------------------------------------------------------
//! \fn void PnetcdfOutput:::WriteOutputFile(Mesh *pm)
//  \brief Cycles over all MeshBlocks and writes OutputData in PNETCDF format
//         One timestep per file

void PnetcdfOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag)
{
  // create filename: "file_basename"+"."+"file_id"+"."+XXXXX+".nc",
  // where XXXXX = 5-digit file_number
  std::string fname;
  char number[6];
  int err;
  sprintf(number,"%05d",output_params.file_number);

  fname.assign(output_params.file_basename);
  fname.append(".");
  fname.append(output_params.file_id);
  fname.append(".");
  fname.append(number);
  fname.append(".nc");

  // 0. reference radius for spherical polar geometry
  Real radius;
  if (COORDINATE_SYSTEM == "spherical_polar") {
    try {
      radius = pin->GetReal("problem", "radius");
    } catch(std::runtime_error& e) {
      std::stringstream msg; 
      msg << "### FATAL ERROR in PnetcdfOutput::WriteOutputFile" << std::endl
          << "Spherical polar coordinate system must define a reference radius in section <problem>" << std::endl
          << "Use radius = XXX " << std::endl;
      ATHENA_ERROR(msg);
    }
  }

  // 1. open file for output
  int ifile;
  err = ncmpi_create(MPI_COMM_WORLD, fname.c_str(), NC_CLOBBER, MPI_INFO_NULL, &ifile);
  ERR

  // 2. coordinate structure
  size_t nx1 = pm->mesh_size.nx1;
  size_t nx2 = pm->mesh_size.nx2;
  size_t nx3 = pm->mesh_size.nx3;
  size_t nx1f = nx1; if (nx1 > 1) nx1f++;
  size_t nx2f = nx2; if (nx2 > 1) nx2f++;
  size_t nx3f = nx3; if (nx3 > 1) nx3f++;

  // 2. define coordinate
  int idt, idx1, idx2, idx3, idx1f, idx2f, idx3f;
  int ndims = 0;

  ncmpi_def_dim(ifile, "time", NC_UNLIMITED, &idt);

  ncmpi_def_dim(ifile, "x1", nx1, &idx1);
  ndims++;
  if (nx1 > 1) {
    ncmpi_def_dim(ifile, "x1f", nx1f, &idx1f);
    ndims++;
  }

  ncmpi_def_dim(ifile, "x2", nx2, &idx2);
  ndims++;
  if (nx2 > 1) {
    ncmpi_def_dim(ifile, "x2f", nx2f, &idx2f);
    ndims++;
  }

  ncmpi_def_dim(ifile, "x3", nx3, &idx3);
  ndims++;
  if (nx3 > 1) {
    ncmpi_def_dim(ifile, "x3f", nx3f, &idx3f);
    ndims++;
  }

  // 3. define variables
  int ivt, ivx1, ivx2, ivx3, ivx1f, ivx2f, ivx3f;

  ncmpi_def_var(ifile, "time", NC_FLOAT, 1, &idt, &ivt);
  ncmpi_put_att_text(ifile, ivt, "axis", 1, "T");
  ncmpi_put_att_text(ifile, ivt, "long_name", 4, "time");
  if (COORDINATE_SYSTEM == "spherical_polar") {
    ncmpi_put_att_text(ifile, ivt, "units", 25, "seconds since 1-1-1 0:0:0");
  } else
    ncmpi_put_att_text(ifile, ivt, "units", 7, "seconds");

  ncmpi_def_var(ifile, "x1", NC_FLOAT, 1, &idx1, &ivx1);
  ncmpi_put_att_text(ifile, ivx1, "axis", 1, "Z");
  if (COORDINATE_SYSTEM == "spherical_polar") {
    ncmpi_put_att_text(ifile, ivx1, "units", 6, "meters");
    ncmpi_put_att_text(ifile, ivx1, "long_name", 8, "altitude");
  } else {
    ncmpi_put_att_text(ifile, ivx1, "units", 6, "meters");
    ncmpi_put_att_text(ifile, ivx1, "long_name", 27, "z-coordinate at cell center");
  }
  if (nx1 > 1) {
    ncmpi_def_var(ifile, "x1f", NC_FLOAT, 1, &idx1f, &ivx1f);
    ncmpi_put_att_text(ifile, ivx1f, "axis", 1, "Z");
    if (COORDINATE_SYSTEM == "spherical_polar") {
      ncmpi_put_att_text(ifile, ivx1f, "units", 6, "meters");
      ncmpi_put_att_text(ifile, ivx1f, "long_name", 8, "altitude");
    } else {
      ncmpi_put_att_text(ifile, ivx1f, "units", 6, "meters");
      ncmpi_put_att_text(ifile, ivx1f, "long_name", 25, "z-coordinate at cell face");
    }
  }

  ncmpi_def_var(ifile, "x2", NC_FLOAT, 1, &idx2, &ivx2);
  ncmpi_put_att_text(ifile, ivx2, "axis", 1, "Y");
  if (COORDINATE_SYSTEM == "spherical_polar") {
    ncmpi_put_att_text(ifile, ivx2, "units", 7, "degrees");
    ncmpi_put_att_text(ifile, ivx2, "long_name", 10, "colatitude");
  } else {
    ncmpi_put_att_text(ifile, ivx2, "units", 6, "meters");
    ncmpi_put_att_text(ifile, ivx2, "long_name", 27, "y-coordinate at cell center");
  }
  if (nx2 > 1) {
    ncmpi_def_var(ifile, "x2f", NC_FLOAT, 1, &idx2f, &ivx2f);
    ncmpi_put_att_text(ifile, ivx2f, "axis", 1, "Y");
    if (COORDINATE_SYSTEM == "spherical_polar") {
      ncmpi_put_att_text(ifile, ivx2f, "units", 7, "degrees");
      ncmpi_put_att_text(ifile, ivx2f, "long_name", 10, "colatitude");
    } else {
      ncmpi_put_att_text(ifile, ivx2f, "units", 6, "meters");
      ncmpi_put_att_text(ifile, ivx2f, "long_name", 25, "y-coordinate at cell face");
    }
  }

  ncmpi_def_var(ifile, "x3", NC_FLOAT, 1, &idx3, &ivx3);
  ncmpi_put_att_text(ifile, ivx3, "axis", 1, "X");
  if (COORDINATE_SYSTEM == "spherical_polar") {
    ncmpi_put_att_text(ifile, ivx3, "units", 7, "degrees");
    ncmpi_put_att_text(ifile, ivx3, "long_name", 9, "longitude");
  } else {
    ncmpi_put_att_text(ifile, ivx3, "units", 6, "meters");
    ncmpi_put_att_text(ifile, ivx3, "long_name", 27, "x-coordinate at cell center");
  }
  if (nx3 > 1) {
    ncmpi_def_var(ifile, "x3f", NC_FLOAT, 1, &idx3f, &ivx3f);
    ncmpi_put_att_text(ifile, ivx3f, "axis", 1, "X");
    if (COORDINATE_SYSTEM == "spherical_polar") {
      ncmpi_put_att_text(ifile, ivx3f, "units", 7, "degrees");
      ncmpi_put_att_text(ifile, ivx3f, "long_name", 9, "longitude");
    } else {
      ncmpi_put_att_text(ifile, ivx3f, "units", 6, "meters");
      ncmpi_put_att_text(ifile, ivx3f, "long_name", 25, "x-coordinate at cell face");
    }
  }

  MeshBlock *pmb = pm->pblock;
  LoadOutputData(pmb);
  OutputData *pdata = pfirst_data_;

  // count total variables (vector variables are expanded into flat scalars)
  int nvars = 0;
  while (pdata != NULL) {
    nvars += pdata->data.GetDim4();
    pdata = pdata->pnext;
  }

  int iax1[2]  = {idt, idx1};
  int iaxis[4]  = {idt, idx3, idx2, idx1};
  int iaxis1[4] = {idt, idx3, idx2, idx1f};
  int iaxis2[4] = {idt, idx3, idx2f, idx1};
  int iaxis3[4] = {idt, idx3f, idx2, idx1};
  int *var_ids = new int [nvars];
  int *ivar = var_ids;

  pdata = pfirst_data_;
  while (pdata != NULL) {
    std::string name;
    for (int n = 0; n < pdata->data.GetDim4(); ++n) {
      size_t pos = pdata->name.find('?');
      if (pdata->data.GetDim4() == 1) { // SCALARS
        if (pos < pdata->name.length()) {  // find '?'
          name = pdata->name.substr(0, pos) + pdata->name.substr(pos + 1);
        } else
          name = pdata->name;
      } else {  // VECTORS
        char c[16]; sprintf(c, "%d", n + 1);
        if (pos < pdata->name.length()) {  // find '?'
          name = pdata->name.substr(0, pos) + c + pdata->name.substr(pos + 1);
        } else
          name = pdata->name + c;
      }
      if (pdata->grid == "CCF")
        ncmpi_def_var(ifile, name.c_str(), NC_FLOAT, 4, iaxis1, ivar);
      else if ((pdata->grid == "CFC") && (nx2 > 1))
        ncmpi_def_var(ifile, name.c_str(), NC_FLOAT, 4, iaxis2, ivar);
      else if ((pdata->grid == "FCC") && (nx3 > 1))
        ncmpi_def_var(ifile, name.c_str(), NC_FLOAT, 4, iaxis3, ivar);
      else if (pdata->grid == "--C")
        ncmpi_def_var(ifile, name.c_str(), NC_FLOAT, 2, iax1, ivar);
      else
        ncmpi_def_var(ifile, name.c_str(), NC_FLOAT, 4, iaxis, ivar);

      if (name == "rho") {
        ncmpi_put_att_text(ifile, *ivar, "units", 6, "kg/m^3");
        ncmpi_put_att_text(ifile, *ivar, "long_name", 7, "density");
      } else if (name == "press") {
        ncmpi_put_att_text(ifile, *ivar, "units", 2, "pa");
        ncmpi_put_att_text(ifile, *ivar, "long_name", 8, "pressure");
      } else if (name == "vel1") {
        ncmpi_put_att_text(ifile, *ivar, "units", 3, "m/s");
        ncmpi_put_att_text(ifile, *ivar, "long_name", 24, "velocity in z-coordinate");
      } else if (name == "vel2") {
        ncmpi_put_att_text(ifile, *ivar, "units", 3, "m/s");
        ncmpi_put_att_text(ifile, *ivar, "long_name", 24, "velocity in y-coordinate");
      } else if (name == "vel3") {
        ncmpi_put_att_text(ifile, *ivar, "units", 3, "m/s");
        ncmpi_put_att_text(ifile, *ivar, "long_name", 24, "velocity in x-coordinate");
      } else if (name == "temp") {
        ncmpi_put_att_text(ifile, *ivar, "units", 1, "K");
        ncmpi_put_att_text(ifile, *ivar, "long_name", 24, "temperature");
      } else if (name == "theta") {
        ncmpi_put_att_text(ifile, *ivar, "units", 1, "K");
        ncmpi_put_att_text(ifile, *ivar, "long_name", 21, "potential temperature");
      } else if (name.substr(name.length()-3) == "tau") {
        ncmpi_put_att_text(ifile, *ivar, "units", 1, "1");
        ncmpi_put_att_text(ifile, *ivar, "long_name", 13, "optical depth");
      } else if (name.substr(name.length()-5) == "flxup") {
        ncmpi_put_att_text(ifile, *ivar, "units", 5, "w/m^2");
        ncmpi_put_att_text(ifile, *ivar, "long_name", 21, "upward radiative flux");
      } else if (name.substr(name.length()-5) == "flxdn") {
        ncmpi_put_att_text(ifile, *ivar, "units", 5, "w/m^2");
        ncmpi_put_att_text(ifile, *ivar, "long_name", 23, "downward radiative flux");
      }

      ivar++;
    }
    pdata = pdata->pnext;
  }

  err = ncmpi_enddef(ifile);
  ERR

  // 4. allocate data buffer, MPI requests and status
  MeshBlock *pb = pmb;
  int max_ncells = 0,  nmb = 0;
  while (pb != NULL) {
    int nf1 = pb->ie - pb->is + 2*NGHOST;
    int nf2 = pb->je - pb->js + 2*NGHOST;
    int nf3 = pb->ke - pb->ks + 2*NGHOST;
    max_ncells  = std::max(max_ncells,  nf1*nf2*nf3);
    nmb++;  // number of meshblocks this rank
    pb = pb->next;
  }
  int nbufs = nmb*(ndims + nvars);
  int *reqs = new int [nbufs];
  int *stts = new int [nbufs];
  float **buf = new float* [nbufs];
  buf[0] = new float [nbufs*max_ncells];
  for (int i = 0; i < nbufs; ++i)
    buf[i] = buf[0] + i*max_ncells;
  int *ir = reqs;
  float **ib = buf;

  // 5. first meshblock writes time
  float time = (float)pm->time;
  MPI_Offset stime = 0, etime = 1;
  err = ncmpi_put_vara_float_all(ifile, ivt, &stime, &etime, &time);
  ERR

  // Loop over MeshBlocks
  //do {
  while (pmb != nullptr) {
    // set start/end array indices depending on whether ghost zones are included
    out_is=pmb->is; out_ie=pmb->ie;
    out_js=pmb->js; out_je=pmb->je;
    out_ks=pmb->ks; out_ke=pmb->ke;
    
    // FIXME: include_ghost zones probably doesn't work with grids other than CCC
    if (output_params.include_ghost_zones) {
      out_is -= NGHOST; out_ie += NGHOST;
      if (out_js != out_je) {out_js -= NGHOST; out_je += NGHOST;}
      if (out_ks != out_ke) {out_ks -= NGHOST; out_ke += NGHOST;}
    }

    // 6. each meshblock writes coordinates and variables
    MPI_Offset ncells1 = out_ie - out_is + 1;
    MPI_Offset ncells2 = out_je - out_js + 1;
    MPI_Offset ncells3 = out_ke - out_ks + 1;
    MPI_Offset nfaces1 = ncells1; if (ncells1 > 1) nfaces1++;
    MPI_Offset nfaces2 = ncells2; if (ncells2 > 1) nfaces2++;
    MPI_Offset nfaces3 = ncells3; if (ncells3 > 1) nfaces3++;

    MPI_Offset start[4] = {0, ncells3*pmb->loc.lx3, ncells2*pmb->loc.lx2, ncells1*pmb->loc.lx1};
    MPI_Offset count[4]  = {1, ncells3, ncells2, ncells1};
    MPI_Offset count1[4] = {1, ncells3, ncells2, nfaces1};
    MPI_Offset count2[4] = {1, ncells3, nfaces2, ncells1};
    MPI_Offset count3[4] = {1, nfaces3, ncells2, ncells1};

    MPI_Offset start_ax1[2] = {0, ncells1*pmb->loc.lx1};
    MPI_Offset count_ax1[2] = {1, ncells1};

    if (COORDINATE_SYSTEM == "spherical_polar")
      for (int i = out_is; i <= out_ie; ++i)
        (*ib)[i-out_is] = (float)(pmb->pcoord->x1v(i)) - radius;
    else
      for (int i = out_is; i <= out_ie; ++i)
        (*ib)[i-out_is] = (float)(pmb->pcoord->x1v(i));
    err = ncmpi_iput_vara_float(ifile, ivx1, start+3, count+3, *ib++, ir++);
    ERR

    if (nx1 > 1) {
      if (COORDINATE_SYSTEM == "spherical_polar")
        for (int i = out_is; i <= out_ie + 1; ++i)
          (*ib)[i-out_is] = (float)(pmb->pcoord->x1f(i)) - radius;
      else
        for (int i = out_is; i <= out_ie + 1; ++i)
          (*ib)[i-out_is] = (float)(pmb->pcoord->x1f(i));
      err = ncmpi_iput_vara_float(ifile, ivx1f, start+3, count1+3, *ib++, ir++);
      ERR
    }

    if (COORDINATE_SYSTEM == "spherical_polar")
      for (int j = out_js; j <= out_je; ++j)
        (*ib)[j-out_js] = (float)rad2deg(pmb->pcoord->x2v(j));
    else
      for (int j = out_js; j <= out_je; ++j)
        (*ib)[j-out_js] = (float)(pmb->pcoord->x2v(j));
    err = ncmpi_iput_vara_float(ifile, ivx2, start+2, count+2, *ib++, ir++);
    ERR

    if (nx2 > 1) {
      if (COORDINATE_SYSTEM == "spherical_polar")
        for (int j = out_js; j <= out_je + 1; ++j)
          (*ib)[j-out_js] = (float)rad2deg(pmb->pcoord->x2f(j));
      else
        for (int j = out_js; j <= out_je + 1; ++j)
          (*ib)[j-out_js] = (float)(pmb->pcoord->x2f(j));
      err = ncmpi_iput_vara_float(ifile, ivx2f, start+2, count2+2, *ib++, ir++);
      ERR
    }

    if (COORDINATE_SYSTEM == "spherical_polar")
      for (int k = out_ks; k <= out_ke; ++k)
        (*ib)[k-out_ks] = (float)rad2deg(pmb->pcoord->x3v(k));
    else
      for (int k = out_ks; k <= out_ke; ++k)
        (*ib)[k-out_ks] = (float)(pmb->pcoord->x3v(k));
    err = ncmpi_iput_vara_float(ifile, ivx3, start+1, count+1, *ib++, ir++);
    ERR

    if (nx3 > 1) {
      if (COORDINATE_SYSTEM == "spherical_polar")
        for (int k = out_ks; k <= out_ke + 1; ++k)
          (*ib)[k-out_ks] = (float)rad2deg(pmb->pcoord->x3f(k));
      else
        for (int k = out_ks; k <= out_ke + 1; ++k)
          (*ib)[k-out_ks] = (float)(pmb->pcoord->x3f(k));
      err = ncmpi_iput_vara_float(ifile, ivx3f, start+1, count3+1, *ib++, ir++);
      ERR
    }

    ivar = var_ids;
    pdata = pfirst_data_;
    while (pdata != NULL) {
      int nvar = pdata->data.GetDim4();
      if (pdata->grid == "CCF") {
        for (int n = 0; n < nvar; n++) {
          float *it = *ib;
          for (int k = out_ks; k <= out_ke; ++k)
            for (int j = out_js; j <= out_je; ++j)
              for (int i = out_is; i <= out_ie+1; ++i)
                *it++ = (float)pdata->data(n,k,j,i);
          err = ncmpi_iput_vara_float(ifile, *ivar++, start, count1, *ib++, ir++);
          ERR
        }
      } else if ((pdata->grid == "CFC") && (nx2 > 1)) {
        for (int n = 0; n < nvar; n++) {
          float *it = *ib;
          for (int k = out_ks; k <= out_ke; ++k)
            for (int j = out_js; j <= out_je+1; ++j)
              for (int i = out_is; i <= out_ie; ++i, ++it)
                *it = (float)pdata->data(n,k,j,i);
          err = ncmpi_iput_vara_float(ifile, *ivar++, start, count2, *ib++, ir++);
          ERR
        }
      } else if ((pdata->grid == "FCC") && (nx3 > 1)) {
        for (int n = 0; n < nvar; n++) {
          float *it = *ib;
          for (int k = out_ks; k <= out_ke+1; ++k)
            for (int j = out_js; j <= out_je; ++j)
              for (int i = out_is; i <= out_ie; ++i, ++it)
                *it = (float)pdata->data(n,k,j,i);
          err = ncmpi_iput_vara_float(ifile, *ivar++, start, count3, *ib++, ir++);
          ERR
        }
      } else if (pdata->grid == "--C") {
        for (int n = 0; n < nvar; n++) {
          float *it = *ib;
          for (int i = out_is; i <= out_ie; ++i, ++it)
            *it = (float)pdata->data(n,i);
          err = ncmpi_iput_vara_float(ifile, *ivar++, start_ax1, count_ax1, *ib++, ir++);
          ERR
        }
      } else {
        for (int n = 0; n < nvar; n++) {
          float *it = *ib;
          for (int k = out_ks; k <= out_ke; ++k)
            for (int j = out_js; j <= out_je; ++j)
              for (int i = out_is; i <= out_ie; ++i, ++it)
                *it = (float)pdata->data(n,k,j,i);
          err = ncmpi_iput_vara_float(ifile, *ivar++, start, count, *ib++, ir++);
          ERR
        }
      }

      pdata = pdata->pnext;
    }

    ClearOutputData();  // required when LoadOutputData() is used.
    pmb = pmb->next;
    if (pmb != nullptr) LoadOutputData(pmb);
  }
  //} while (LoadOutputData(pmb));   // end loop over MeshBLocks

  // 7. wait for all writings to complete
  ncmpi_wait_all(ifile, nbufs, reqs, stts);
  for (int i = 0; i < nbufs; ++i)
    if (stts[i] != NC_NOERR) {
      std::stringstream msg; 
      msg << "### FATAL ERROR in PnetcdfOutput::WriteOutputFile"
          << std::endl << "Nonblocking write fails on request "
          << i << std::endl
          << ncmpi_strerror(stts[i]);
      throw std::runtime_error(msg.str().c_str());
    }

  // 8. close nc file
  ncmpi_close(ifile);

  // 9. clear data buffers
  delete [] buf[0];
  delete [] buf;
  delete [] var_ids;
  delete [] reqs;
  delete [] stts;

  // 10. increment counters
  output_params.file_number++;
  output_params.next_time += output_params.dt;
  pin->SetInteger(output_params.block_name, "file_number", output_params.file_number);
  pin->SetReal(output_params.block_name, "next_time", output_params.next_time);

  return;
}

#endif // PNETCDFOUTPUT
