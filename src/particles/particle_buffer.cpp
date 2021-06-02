// C/C++ headers
#include <sstream>
#include <cstddef>
#include <functional>
#include <iostream>
#ifdef MPI_PARALLEL
  // defined in particles.cpp
  extern MPI_Datatype MPI_PARTICLE;
#endif

// Athena++ classes headers
#include "../mesh/mesh.hpp"
#include "../globals.hpp"
#include "../bvals/bvals.hpp"
#include "material_point.hpp"
#include "particle_buffer.hpp"
#include "particles.hpp"

ParticleBuffer::ParticleBuffer(Particles *ppar):
  pmy_particle(ppar)
{
  for (int i = 0; i < 56; ++i) {
    particle_flag_[i] = BoundaryStatus::completed;

#ifdef MPI_PARALLEL
    req_particle_send_[i]=MPI_REQUEST_NULL;
    req_particle_recv_[i]=MPI_REQUEST_NULL;
#endif
  }

  MeshBlock *pmb = ppar->pmy_block;
  for (int n = 0; n < pmb->pbval->nneighbor; ++n) {
    NeighborBlock &nb = pmb->pbval->neighbor[n];
    particle_flag_[nb.bufid] = BoundaryStatus::waiting;
  }
}

ParticleBuffer::~ParticleBuffer() {}

int ParticleBuffer::CreateMPITag(int lid, int tid) {
  int TAG_PARTICLE = 15;
  int tag = BoundaryBase::CreateBvalsMPITag(lid, tid, TAG_PARTICLE);
  std::string str = pmy_particle->myname + std::to_string(tag);
  return std::hash<std::string>{}(str)%(1<<24);
}

void ParticleBuffer::DetachParticle(std::vector<MaterialPoint> &mp)
{
  for (size_t i = 0; i < bufid.size(); ++i)
    if (bufid[i] >= 0) { // if particle is still alive
      size_t j = mp.size() - i - 1;
      particle_send_[bufid[i]].push_back(mp[j]);
    }

  mp.resize(mp.size() - bufid.size());
  bufid.clear();
}

void ParticleBuffer::SendParticle()
{
  for (int n = 0; n < pmy_particle->pmy_block->pbval->nneighbor; ++n) {
    NeighborBlock &nb = pmy_particle->pmy_block->pbval->neighbor[n];

    if (nb.snb.rank == Globals::my_rank) {  // on the same process
      MeshBlock *pmb = pmy_particle->pmy_block->pmy_mesh->FindMeshBlock(nb.snb.gid);
      pmb->ppar->ppb->particle_recv_[nb.targetid] = particle_send_[nb.bufid];
      pmb->ppar->ppb->particle_flag_[nb.targetid] = BoundaryStatus::arrived;
    }
#ifdef MPI_PARALLEL
    else { // MPI
      int tag = CreateBvalsMPITag(nb.lid, nb.targetid);
      int ssize = particle_send_[nb.bufid].size();
      MPI_Isend(particle_send_[nb.bufid].data(), ssize, MPI_PARTICLE,
                nb.snb.rank, tag, MPI_COMM_WORLD, &req_particle_send_[nb.bufid]);
    }
#endif
  }
}

void ParticleBuffer::RecvParticle()
{
  int rsize, tag;

  MeshBlock *pmb = pmy_particle->pmy_block;
#ifdef MPI_PARALLEL
  MPI_Status status;
  for (int n = 0; n < pmb->pbval->nneighbor; ++n) {
    NeighborBlock &nb = pmb->pbval->neighbor[n];
    if (nb.snb.rank == Globals::my_rank) continue; // local boundary received

    if (particle_flag_[nb.bufid] == BoundaryStatus::waiting) {
      int tag = CreateBvalsMPITag(pmb->lid, nb.bufid);
      MPI_Probe(nb.rank, tag, MPI_COMM_WORLD, &status);
      MPI_Get_count(&status, MPI_PARTICLE, &rsize);
      particle_recv_[nb.bufid].resize(rsize);
      MPI_Irecv(particle_recv_[nb.bufid].data(), rsize, MPI_PARTICLE,
                nb.rank, tag, MPI_COMM_WORLD, &req_particle_recv_[nb.bufid]);
    }
  }
#endif
}

// attach particle to the particle chain
bool ParticleBuffer::AttachParticle(std::vector<MaterialPoint>& mp)
{
  bool success = true;
  int test;

  MeshBlock *pmb = pmy_particle->pmy_block;
  Mesh *pm = pmb->pmy_mesh;
  for (int n = 0; n < pmb->pbval->nneighbor; ++n) {
    NeighborBlock &nb = pmb->pbval->neighbor[n];
    if (particle_flag_[nb.bufid] == BoundaryStatus::completed) continue;

#ifdef MPI_PARALLEL
    if (particle_flag_[nb.bufid] == BoundaryStatus::waiting) {
      MPI_Test(&req_particle_recv_[nb.bufid], &test, MPI_STATUS_IGNORE);
      if (test)
        particle_flag_[nb.bufid] = BoundaryStatus::arrived;
      else {
        success = false;
        continue;
      }
    }
#endif

    if (particle_flag_[nb.bufid] == BoundaryStatus::arrived) {
      std::vector<MaterialPoint>::iterator it = particle_recv_[nb.bufid].begin();
      for (; it != particle_recv_[nb.bufid].end(); ++it) {
        // 0:INNER_X1, 1:OUTER_X1
        if (pm->mesh_bcs[nb.ni.ox1+1>>1] == BoundaryFlag::periodic)
          it->x1 += nb.ni.ox1*(pm->mesh_size.x1max - pm->mesh_size.x1min);
        // 2:INNER_X2, 3:OUTER_X2
        if (pm->mesh_bcs[2+(nb.ni.ox2+1>>1)] == BoundaryFlag::periodic)
          it->x2 += nb.ni.ox2*(pm->mesh_size.x2max - pm->mesh_size.x2min);
        // 4:INNER_X3, 5:OUTER_X3
        if (pm->mesh_bcs[4+(nb.ni.ox3+1>>1)] == BoundaryFlag::periodic) 
          it->x3 += nb.ni.ox3*(pm->mesh_size.x3max - pm->mesh_size.x3min);

        bool out_of_domain = (it->x1 < pmb->block_size.x1min || it->x1 > pmb->block_size.x1max) ||
             (pm->f2 && (it->x2 < pmb->block_size.x2min || it->x2 > pmb->block_size.x2max)) ||
             (pm->f3 && (it->x3 < pmb->block_size.x3min || it->x3 > pmb->block_size.x3max));
        if (out_of_domain) {
          success = false;
          std::stringstream msg;
          msg << "### FATAL ERROR in ParticleBuffer::AttachParticles. Particles " 
              << pmy_particle->myname
              << " moved out of MeshBlock limits" << std::endl;
          msg << it->x1 << " " << pmb->block_size.x1min << " " << pmb->block_size.x1max << std::endl;
          msg << it->x2 << " " << pmb->block_size.x2min << " " << pmb->block_size.x2max << std::endl;
          msg << it->x3 << " " << pmb->block_size.x3min << " " << pmb->block_size.x3max << std::endl;
          ATHENA_ERROR(msg);
        }

        mp.push_back(*it);
      }
      particle_flag_[nb.bufid] = BoundaryStatus::completed; // completed
    }
  }

  return success;
}

void ParticleBuffer::ClearBoundary()
{
  MeshBlock *pmb = pmy_particle->pmy_block;
  for (int n = 0; n < pmb->pbval->nneighbor; ++n) {
    NeighborBlock &nb = pmb->pbval->neighbor[n];
    particle_flag_[nb.bufid] = BoundaryStatus::waiting;
    particle_send_[nb.bufid].clear();
    particle_recv_[nb.bufid].clear();

#ifdef MPI_PARALLEL
    if (nb.snb.rank != Globals::my_rank)
      MPI_Wait(&req_particle_send_[nb.bufid], MPI_STATUS_IGNORE);
#endif
  }
}
