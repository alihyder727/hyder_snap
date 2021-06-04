/** @file particles_table.hpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Thursday June 3, 2021 14:23:40 PDT
 * @bug No known bugs.
 */

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
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "outputs.hpp"
#include "../particles/material_point.hpp"
#include "../particles/particles.hpp"

//----------------------------------------------------------------------------------------
// ParticlesTableOutput constructor
// destructor not required for this derived class

ParticlesTableOutput::ParticlesTableOutput(OutputParameters oparams)
  : OutputType(oparams)
{}

//----------------------------------------------------------------------------------------
//! \fn void ParticlesTableOutput:::WriteOutputFile(Mesh *pm)
//  \brief writes OutputData to file in tabular format using C style fprintf
//         Writes one file per ParticleGroup per MeshBlock

void ParticlesTableOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag)
{
  MeshBlock *pmb=pm->pblock;

  // Loop over MeshBlocks
  while (pmb != nullptr) {
    Particles *ppar = pmb->ppar;
    // Loop over ParticleGroup
    while (ppar != nullptr) {
      // create filename: "file_basename"+"."+"name"+"."+"blockid"+"."+?????+".ptab",
      // where ????? = 5-digit file_number
      std::string fname;
      char number[6];
      sprintf(number,"%05d",output_params.file_number);
      char blockid[12];
      sprintf(blockid,"block%d",pmb->gid);
  
      fname.assign(output_params.file_basename);
      fname.append(".");
      fname.append(ppar->myname);
      fname.append(".");
      fname.append(blockid);
      fname.append(".");
      fname.append(number);
      fname.append(".ptab");

      // open file for output
      FILE *pfile;
      std::stringstream msg;
      if ((pfile = fopen(fname.c_str(),"w")) == nullptr){
        msg << "### FATAL ERROR in function [ParticlesTableOutput::WriteOutputFile]"
            <<std::endl<< "Output file '" <<fname<< "' could not be opened" <<std::endl;
        throw std::runtime_error(msg.str().c_str());
      }

      // print file header
      fprintf(pfile,"# Athena++ data at time=%e",pm->time);
      fprintf(pfile,"  cycle=%d",pmb->pmy_mesh->ncycle);
      fprintf(pfile,"  particle=%s",ppar->myname.c_str());
      fprintf(pfile,"  number of particles=%ld \n",ppar->mp.size());

      // write x1, x2, x3 column headers
      fprintf(pfile,"#");
      fprintf(pfile,"%5s","cat");
      fprintf(pfile,"%16s","mass");
      fprintf(pfile,"%16s","x1");
      fprintf(pfile,"%16s","x2");
      fprintf(pfile,"%16s","x3");
      fprintf(pfile,"%16s","v1");
      fprintf(pfile,"%16s","v2");
      fprintf(pfile,"%16s","v3");
      for (int j = 0; j < NREAL_PARTICLE_DATA; ++j)
        fprintf(pfile, "%14s%02d", "RDATA", j);
      for (int j = 0; j < NINT_PARTICLE_DATA; ++j)
        fprintf(pfile, "%14s%02d", "IDATA", j);
      fprintf(pfile,"\n"); // terminate line

      // loop over all particles
      std::vector<MaterialPoint>::iterator it = ppar->mp.begin();
      for (; it != ppar->mp.end(); ++it) {
        fprintf(pfile, "%-6d", it->ct);
        fprintf(pfile, output_params.data_format.c_str(), it->mass);
        fprintf(pfile, output_params.data_format.c_str(), it->x1);
        fprintf(pfile, output_params.data_format.c_str(), it->x2);
        fprintf(pfile, output_params.data_format.c_str(), it->x3);
        fprintf(pfile, output_params.data_format.c_str(), it->v1);
        fprintf(pfile, output_params.data_format.c_str(), it->v2);
        fprintf(pfile, output_params.data_format.c_str(), it->v3);
        #if NREAL_PARTICLE_DATA > 0
          for (size_t j = 0; j < NREAL_PARTICLE_DATA; ++j)
            fprintf(pfile, output_params.data_format.c_str(), it->rr[j]);
        #endif

        #if NINT_PARTICLE_DATA > 0
          for (size_t j = 0; j < NINT_PARTICLE_DATA; ++j)
            fprintf(pfile, output_params.data_format.c_str(), it->ii[j]);
        #endif
        fprintf(pfile,"\n"); // terminate line
      }

      // close file, and next variable
      fclose(pfile);
      ppar = ppar->next;
    } // end loop over ParticleGroup
    pmb = pmb->next;
  }  // end loop over MeshBlocks

  // increment counters
  output_params.file_number++;
  output_params.next_time += output_params.dt;
  pin->SetInteger(output_params.block_name, "file_number", output_params.file_number);
  pin->SetReal(output_params.block_name, "next_time", output_params.next_time);

  return;
}
