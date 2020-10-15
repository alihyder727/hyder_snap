//! \file netcdf.cpp
//  \brief writes output data in NETCDF format.
//  Data is written in RECTILINEAR_GRID geometry, in BINARY format, and in FLOAT type
//  Writes one file per MeshBlock.

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
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../coordinates/coordinates.hpp"
#include "outputs.hpp"

// Only proceed if NETCDF output enabled
#ifdef NETCDFOUTPUT

// External library headers
#include <netcdf.h> // nc_[create|close], nc_enddef, nc_strerror,
                    // nc_def_[dim|var], nc_put_var_float, nc_put_vara_float

//----------------------------------------------------------------------------------------
// NetcdfOutput constructor
// destructor - not needed for this derived class

NetcdfOutput::NetcdfOutput(OutputParameters oparams)
  : OutputType(oparams)
{
}

//----------------------------------------------------------------------------------------
//! \fn void NetcdfOutput:::WriteOutputFile(Mesh *pm)
//  \brief Cycles over all MeshBlocks and writes OutputData in NETCDF format, one
//         MeshBlock per file

void NetcdfOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag)
{
  // find out total number of blocks
//  int nbtotal;
//#ifdef MPI_PARALLEL
//  MPI_Allreduce(&pm->nbtotal, &nbtotal, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
//#else
//  nbtotal = pm->nbtotal;
//#endif

  MeshBlock *pmb=pm->pblock;

  // Loop over MeshBlocks
  while (pmb != NULL) {
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

    // set ptrs to data in OutputData linked list, then slice/sum as needed
    LoadOutputData(pmb);
    if (TransformOutputData(pmb) == false) {
      ClearOutputData();  // required when LoadOutputData() is used.
      pmb=pmb->next;
      continue;
    } // skip if slice was out of range

    // create filename: "file_basename"+ "."+"blockid"+"."+"file_id"+"."+XXXXX+".nc",
    // where XXXXX = 5-digit file_number
    std::string fname;
    char number[6];
    sprintf(number,"%05d",output_params.file_number);
    char blockid[12];
    sprintf(blockid,"block%d",pmb->gid);

    fname.assign(output_params.file_basename);
    fname.append(".");
    fname.append(blockid);
    fname.append(".");
    fname.append(output_params.file_id);
    fname.append(".");
    fname.append(number);
    fname.append(".nc");

    // 1. open file for output
    std::stringstream msg;
    int ifile;

    nc_create(fname.c_str(), NC_CLASSIC_MODEL, &ifile);

    // 2. coordinate structure
    int ncells1 = out_ie - out_is + 1;
    int ncells2 = out_je - out_js + 1;
    int ncells3 = out_ke - out_ks + 1;
    int nfaces1 = ncells1; if (ncells1 > 1) nfaces1++;
    int nfaces2 = ncells2; if (ncells2 > 1) nfaces2++;
    int nfaces3 = ncells3; if (ncells3 > 1) nfaces3++;

    // 2. define coordinate
    int idt, idx1, idx2, idx3, idx1f, idx2f, idx3f;
    // time
    nc_def_dim(ifile, "time", NC_UNLIMITED, &idt);

    nc_def_dim(ifile, "x1", ncells1, &idx1);
    if (ncells1 > 1)
      nc_def_dim(ifile, "x1f", nfaces1, &idx1f);

    nc_def_dim(ifile, "x2", ncells2, &idx2);
    if (ncells2 > 1)
      nc_def_dim(ifile, "x2f", nfaces2, &idx2f);

    nc_def_dim(ifile, "x3", ncells3, &idx3);
    if (ncells3 > 1)
      nc_def_dim(ifile, "x3f", nfaces3, &idx3f);

    // 3. define variables
    int ivt, ivx1, ivx2, ivx3, ivx1f, ivx2f, ivx3f;
    int loc[4] = {(int)pmb->loc.lx1, (int)pmb->loc.lx2, (int)pmb->loc.lx3, pmb->loc.level};
    int pos[4];

    nc_def_var(ifile, "time", NC_FLOAT, 1, &idt, &ivt);
    nc_put_att_text(ifile, ivt, "axis", 1, "T");
    nc_put_att_text(ifile, ivt, "units", 1, "s");

    nc_def_var(ifile, "x1", NC_FLOAT, 1, &idx1, &ivx1);
    nc_put_att_text(ifile, ivx1, "axis", 1, "Z");
    nc_put_att_text(ifile, ivx1, "units", 1, "m");

    pos[0] = 1; pos[1] = pmb->pmy_mesh->mesh_size.nx1;
    pos[2] = ncells1*loc[0] + 1; pos[3] = ncells1*(loc[0] + 1);
    nc_put_att_int(ifile, ivx1, "domain_decomposition", NC_INT, 4, pos);
    if (ncells1 > 1) {
      nc_def_var(ifile, "x1f", NC_FLOAT, 1, &idx1f, &ivx1f);
      pos[0]--; pos[2]--;
      nc_put_att_int(ifile, ivx1f, "domain_decomposition", NC_INT, 4, pos);
    }

    nc_def_var(ifile, "x2", NC_FLOAT, 1, &idx2, &ivx2);
    nc_put_att_text(ifile, ivx2, "axis", 1, "Y");
    nc_put_att_text(ifile, ivx2, "units", 1, "m");

    pos[0] = 1; pos[1] = pmb->pmy_mesh->mesh_size.nx2;
    pos[2] = ncells2*loc[1] + 1; pos[3] = ncells2*(loc[1] + 1);
    nc_put_att_int(ifile, ivx2, "domain_decomposition", NC_INT, 4, pos);
    if (ncells2 > 1) {
      nc_def_var(ifile, "x2f", NC_FLOAT, 1, &idx2f, &ivx2f);
      pos[0]--; pos[2]--;
      nc_put_att_int(ifile, ivx2f, "domain_decomposition", NC_INT, 4, pos);
    }

    nc_def_var(ifile, "x3", NC_FLOAT, 1, &idx3, &ivx3);
    nc_put_att_text(ifile, ivx3, "axis", 1, "X");
    nc_put_att_text(ifile, ivx3, "units", 1, "m");

    pos[0] = 1; pos[1] = pmb->pmy_mesh->mesh_size.nx3;
    pos[2] = ncells3*loc[2] + 1; pos[3] = ncells3*(loc[2] + 1);
    nc_put_att_int(ifile, ivx3, "domain_decomposition", NC_INT, 4, pos);
    if (ncells3 > 1) {
      nc_def_var(ifile, "x3f", NC_FLOAT, 1, &idx3f, &ivx3f);
      pos[0]--; pos[2]--;
      nc_put_att_int(ifile, ivx3f, "domain_decomposition", NC_INT, 4, pos);
    }

    nc_put_att_int(ifile, NC_GLOBAL, "NumFilesInSet", NC_INT, 1, &pm->nbtotal);
    if ((output_params.variable.compare("diag") == 0) ||
        (output_params.variable.compare("flux") == 0)) {
      float dt = (float)pm->dt/2;
      nc_put_att_float(ifile, NC_GLOBAL, "Time_shift_for_DiagFlux_variables", NC_FLOAT, 1, &dt);
    }

    OutputData *pdata = pfirst_data_;

    // count total variables (vector variables are expanded into flat scalars)
    int total_vars = 0;
    while (pdata != NULL) {
      if (pdata->grid == "-CC")
        total_vars += pdata->data.GetDim3();
      else
        total_vars += pdata->data.GetDim4();
      pdata = pdata->pnext;
    }

    int iax1[2]  = {idt, idx1};
    int iaxis[4]  = {idt, idx3 , idx2 , idx1};
    int iaxis1[4] = {idt, idx3 , idx2 , idx1f};
    int iaxis2[4] = {idt, idx3 , idx2f, idx1};
    int iaxis3[4] = {idt, idx3f, idx2 , idx1};
    int iaxis4[4] = {idt, idx3, idx2};
    int *var_ids = new int [total_vars];
    int *ivar = var_ids;

    pdata = pfirst_data_;
    while (pdata != NULL) {
      std::string name;
      int nvar;
      if (pdata->grid == "-CC") // FIXME: Fix this in radiation output
        nvar = pdata->data.GetDim3();
      else
        nvar = pdata->data.GetDim4();
      for (int n = 0; n < nvar; ++n) {
        if (nvar == 1) { // SCALARS
          size_t pos = pdata->name.find('?');
          if (pos < pdata->name.length()) {  // find '?'
            name = pdata->name.substr(0, pos) + pdata->name.substr(pos + 1);
          } else
            name = pdata->name;
        } else {  // VECTORS
          char c[16]; sprintf(c, "%d", n + 1);
          size_t pos = pdata->name.find('?');
          if (pos < pdata->name.length()) {  // find '?'
            name = pdata->name.substr(0, pos) + c + pdata->name.substr(pos + 1);
          } else
            name = pdata->name + c;
        }
        if (pdata->grid == "CCF")
          nc_def_var(ifile, name.c_str(), NC_FLOAT, 4, iaxis1, ivar++);
        else if ((pdata->grid == "CFC") && (ncells2 > 1))
          nc_def_var(ifile, name.c_str(), NC_FLOAT, 4, iaxis2, ivar++);
        else if ((pdata->grid == "FCC") && (ncells3 > 1))
          nc_def_var(ifile, name.c_str(), NC_FLOAT, 4, iaxis3, ivar++);
        else if (pdata->grid == "--C")
          nc_def_var(ifile, name.c_str(), NC_FLOAT, 2, iax1, ivar++);
        else if (pdata->grid == "-CC")
          nc_def_var(ifile, name.c_str(), NC_FLOAT, 3, iaxis4, ivar++);
        else
          nc_def_var(ifile, name.c_str(), NC_FLOAT, 4, iaxis, ivar++);
      }
      pdata = pdata->pnext;
    }

    nc_enddef(ifile);

    // 4. write variables
    float *data = new float[nfaces1*nfaces2*nfaces3];
    size_t start[4] = {0, 0, 0, 0};
    size_t count[4]  = {1, (size_t)ncells3, (size_t)ncells2, (size_t)ncells1};
    size_t count1[4] = {1, (size_t)ncells3, (size_t)ncells2, (size_t)nfaces1};
    size_t count2[4] = {1, (size_t)ncells3, (size_t)nfaces2, (size_t)ncells1};
    size_t count3[4] = {1, (size_t)nfaces3, (size_t)ncells2, (size_t)ncells1};
    size_t count4[3] = {1, (size_t)ncells3, (size_t)ncells2};
    size_t count_ax1[2] = {1, (size_t)ncells1};

    float time = (float)pm->time;
    nc_put_vara_float(ifile, ivt, start, count, &time);

    for (int i = out_is; i <= out_ie; ++i)
      data[i-out_is] = (float)(pmb->pcoord->x1v(i));
    nc_put_var_float(ifile, ivx1, data);

    if (ncells1 > 1) {
      for (int i = out_is; i <= out_ie + 1; ++i)
        data[i-out_is] = (float)(pmb->pcoord->x1f(i));
      nc_put_var_float(ifile, ivx1f, data);
    }

    for (int j = out_js; j <= out_je; ++j)
      data[j-out_js] = (float)(pmb->pcoord->x2v(j));
    nc_put_var_float(ifile, ivx2, data);

    if (ncells2 > 1) {
      for (int j = out_js; j <= out_je + 1; ++j)
        data[j-out_js] = (float)(pmb->pcoord->x2f(j));
      nc_put_var_float(ifile, ivx2f, data);
    }

    for (int k = out_ks; k <= out_ke; ++k)
      data[k-out_ks] = (float)(pmb->pcoord->x3v(k));
    nc_put_var_float(ifile, ivx3, data);

    if (ncells3 > 1) {
      for (int k = out_ks; k <= out_ke + 1; ++k)
        data[k-out_ks] = (float)(pmb->pcoord->x1f(k));
      nc_put_var_float(ifile, ivx3f, data);
    }

    ivar = var_ids;
    pdata = pfirst_data_;
    while (pdata != NULL) {
      int nvar;
      if (pdata->grid == "-CC")
        nvar = pdata->data.GetDim3();
      else
        nvar = pdata->data.GetDim4();
      if (pdata->grid == "CCF") {
        for (int n = 0; n < nvar; n++) {
          float *it = data;
          for (int k = out_ks; k <= out_ke; ++k)
            for (int j = out_js; j <= out_je; ++j)
              for (int i = out_is; i <= out_ie+1; ++i, ++it)
                *it = (float)pdata->data(n,k,j,i);
          nc_put_vara_float(ifile, *(ivar++), start, count1, data);
        }
      } else if ((pdata->grid == "CFC") && (ncells2 > 1)) {
        for (int n = 0; n < nvar; n++) {
          float *it = data;
          for (int k = out_ks; k <= out_ke; ++k)
            for (int j = out_js; j <= out_je+1; ++j)
              for (int i = out_is; i <= out_ie; ++i, ++it)
                *it = (float)pdata->data(n,k,j,i);
          nc_put_vara_float(ifile, *(ivar++), start, count2, data);
        }
      } else if ((pdata->grid == "FCC") && (ncells3 > 1)) {
        for (int n = 0; n < nvar; n++) {
          float *it = data;
          for (int k = out_ks; k <= out_ke+1; ++k)
            for (int j = out_js; j <= out_je; ++j)
              for (int i = out_is; i <= out_ie; ++i, ++it)
                *it = (float)pdata->data(n,k,j,i);
          nc_put_vara_float(ifile, *(ivar++), start, count3, data);
        }
      } else if (pdata->grid == "--C") {
        for (int n = 0; n < nvar; n++) {
          float *it = data;
          for (int i = out_is; i <= out_ie; ++i, ++it)
            *it = (float)pdata->data(n,i);
          nc_put_vara_float(ifile, *(ivar++), start, count_ax1, data);
        }
      } else if (pdata->grid == "-CC") {
        for (int n = 0; n < nvar; n++) {
          float *it = data;
          for (int k = out_ks; k <= out_ke; ++k)
            for (int j = out_js; j <= out_je; ++j, ++it)
              *it = (float)pdata->data(n,k,j);
          nc_put_vara_float(ifile, *(ivar++), start, count4, data);
        }
      } else {
        for (int n = 0; n < nvar; n++) {
          float *it = data;
          for (int k = out_ks; k <= out_ke; ++k)
            for (int j = out_js; j <= out_je; ++j)
              for (int i = out_is; i <= out_ie; ++i, ++it)
                *it = (float)pdata->data(n,k,j,i);
          nc_put_vara_float(ifile, *(ivar++), start, count, data);
        }
      }

      // doesn't work
      //nc_put_att_text(ifile, *(ivar-1), "output",
      //  output_params.variable.length(), output_params.variable.c_str());
      pdata = pdata->pnext;
    }

    // 5. close nc file
    nc_close(ifile);

    ClearOutputData();  // required when LoadOutputData() is used.
    delete [] data;
    delete [] var_ids;

    pmb=pmb->next;
  }  // end loop over MeshBlocks

  // increment counters
  output_params.file_number++;
  output_params.next_time += output_params.dt;
  pin->SetInteger(output_params.block_name, "file_number", output_params.file_number);
  pin->SetReal(output_params.block_name, "next_time", output_params.next_time);

  return;
}

#endif // NETCDFOUTPUT
