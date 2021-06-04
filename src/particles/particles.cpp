/** @file particles.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Sunday May 30, 2021 13:10:12 PDT
 * @bug No known bugs.
 */

// C++ headers
#include <sstream>
#include <stdexcept>
#ifdef MPI_PARALLEL
#include <mpi.h>
  MPI_Datatype MPI_PARTICLE;
#endif

// Athena++ headers
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "../math/interpolation.h" // interpn, locate
#include "material_point.hpp"
#include "particles.hpp"
#include "particle_buffer.hpp"

Particles::Particles(MeshBlock *pmb, ParameterInput *pin):
  pmy_block(pmb), myname("HEAD"), prev(nullptr), next(nullptr)
{
  ppb = new ParticleBuffer(this);

  char particle_names[1024], *p;
  std::string str = pin->GetOrAddString("particles", "particles", "");
  std::strcpy(particle_names, str.c_str());
  p = std::strtok(particle_names, " ,");

  while (p != NULL) {
    std::stringstream msg;
    std::string name;
    char *c = std::strchr(p, '.');
    if (c != NULL) name = c+1;
    else name = p;
    if (std::strncmp(p, "2pcp", 4) == 0) {
      AddParticles(TwoPhaseCloudParticles(pmb, pin, name));
    } else {
      msg << "### FATAL ERROR in function Particles::Particles"
          << std::endl << "Particles '" << p << "' "
          << "unrecognized." << std::endl;
      ATHENA_ERROR(msg);
    }
    p = std::strtok(NULL, " ,");
  }

#ifdef MPI_PARALLEL
  // MPI_PARTICLE
  int counts[2] = {5+NINT_PARTICLE_DATA, 7+NREAL_PARTICLE_DATA};
  MPI_Datatype types[2] = {MPI_INT, MPI_ATHENA_REAL};
  MPI_Aint disps[2] = {offsetof(MaterialPoint, alive), offsetof(MaterialPoint, mass)};

  MPI_Type_create_struct(2, counts, disps, types, &MPI_PARTICLE);
  MPI_Type_commit(&MPI_PARTICLE);
#endif
}

// constructor, initializes data structure and parameters
Particles::Particles(MeshBlock *pmb, ParameterInput *pin, std::string name):
  pmy_block(pmb), myname(name), prev(nullptr), next(nullptr)
{
  int nc1 = pmb->ncells1, nc2 = pmb->ncells2, nc3 = pmb->ncells3;

  coordinates_.resize(nc3+nc2+nc1);
  for (int k = 0; k < nc3; ++k)
    coordinates_[k] = pmb->pcoord->x3v(k);
  for (int j = 0; j < nc2; ++j)
    coordinates_[nc3+j] = pmb->pcoord->x2v(j);
  for (int i = 0; i < nc1; ++i)
    coordinates_[nc3+nc2+i] = pmb->pcoord->x1v(i);

  lengths_.resize(3);
  lengths_[0] = nc3;
  lengths_[1] = nc2;
  lengths_[2] = nc1;
  vol_.NewAthenaArray(nc1);

  ppb = new ParticleBuffer(this);
}

// destructor
Particles::~Particles()
{
  if (prev != nullptr) prev->next = next;
  if (next != nullptr) next->prev = prev;
  delete ppb;
#ifdef MPI_PARALLEL
  if ((prev == nullptr) && (next == nullptr))
    MPI_Type_free(&MPI_PARTICLE);
#endif
}

Particles::Particles(Particles const& other):
  c(other.c), dc(other.dc), vol_(other.vol_)
{
  if (this == &other) return;
  *this = other;
}

Particles& Particles::operator=(Particles const& other)
{
  pmy_block = other.pmy_block;
  myname = other.myname;
  prev = other.prev;
  next = other.next;
  // assignment operator of AthenaArray does not allocate memory
  c = other.c;
  dc = other.dc;
  mp = other.mp;
  mp1 = other.mp1;

  vol_ = other.vol_;
  coordinates_ = other.coordinates_;
  lengths_ = other.lengths_;
  cnames_ = other.cnames_;

  ppb = new ParticleBuffer(this);

  return *this;
}

// functions
Particles* Particles::FindParticle(std::string name)
{
  std::stringstream msg;
  Particles *p = this;

  while ((p != nullptr) && (p->myname != name)) p = p->next;
  if (p == nullptr) {
    msg << "### FATAL ERROR in Particles::FindParticles"
        << "Particles " << name << " not found" << std::endl;
    ATHENA_ERROR(msg);
  }

  return p;
}

void Particles::ExchangeHydro(AthenaArray<Real> &du, AthenaArray<Real> const &w)
{
  AthenaArray<Real> v1, v2, v3;
  Real loc[3], vel;

  v1.InitWithShallowSlice(const_cast<AthenaArray<Real>&>(w),4,IM1,1);
  v2.InitWithShallowSlice(const_cast<AthenaArray<Real>&>(w),4,IM2,1);
  v3.InitWithShallowSlice(const_cast<AthenaArray<Real>&>(w),4,IM3,1);

  for (std::vector<MaterialPoint>::iterator q = mp.begin(); q != mp.end(); ++q) {
    loc[0] = q->x3;
    loc[1] = q->x2;
    loc[2] = q->x1;

    interpn(&vel, loc, v1.data(), coordinates_.data(), lengths_.data(), 3);
    q->v1 = vel;

    if (lengths_[1] > 1) {
      interpn(&vel, loc, v2.data(), coordinates_.data(), lengths_.data(), 3);
      q->v2 = vel;
    }

    if (lengths_[0] > 1) {
      interpn(&vel, loc, v3.data(), coordinates_.data(), lengths_.data(), 3);
      q->v3 = vel;
    }
  }
}

/*void Particles::SetMeshLocation()
{
  MeshBlock *pmb = pmy_block;
  Real loc[3];
  for (std::vector<MaterialPoint>::iterator q = mp.begin(); q != mp.end(); ++q) {
    loc[0] = q->x3;
    loc[1] = q->x2;
    loc[2] = q->x1;

    if (lengths_[0] > 1)
      q->ck = pmb->ks + locate(coordinates_.data(), loc[0], lengths_[0]);
    if (lengths_[1] > 1)
      q->cj = pmb->js + locate(coordinates_.data() + lengths_[0], loc[1], lengths_[1]);
    q->ci = pmb->is + locate(coordinates_.data() + lengths_[0] + lengths_[1], 
      loc[2], lengths_[2]);
  }
}*/

void Particles::AggregateMass(AthenaArray<Real> &c_sum)
{
  MeshBlock *pmb = pmy_block;
  Real loc[3];
  c_sum.ZeroClear();
  for (std::vector<MaterialPoint>::iterator q = mp.begin(); q != mp.end(); ++q) {
    loc[0] = q->x3;
    loc[1] = q->x2;
    loc[2] = q->x1;

    if (lengths_[0] > 1)
      q->ck = pmb->ks + locate(coordinates_.data(), loc[0], lengths_[0]);
    if (lengths_[1] > 1)
      q->cj = pmb->js + locate(coordinates_.data() + lengths_[0], loc[1], lengths_[1]);
    q->ci = pmb->is + locate(coordinates_.data() + lengths_[0] + lengths_[1], 
      loc[2], lengths_[2]);
    c_sum(q->ct, q->ck, q->cj, q->ci) += q->mass;
  }

  for (int t = 0; t < c_sum.GetDim4(); ++t)
    for (int k = pmb->ks; k <= pmb->ke; ++k)
      for (int j = pmb->js; j <= pmb->je; ++j) {
        pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol_);
        for (int i = pmb->is; i <= pmb->ie; ++i)
          c_sum(t,k,j,i) /= vol_(i);
      }
}

void Particles::Particulate(AthenaArray<Real> &c_dif)
{
  MeshBlock *pmb = pmy_block;
  Coordinates *pco = pmb->pcoord;
  for (int t = 0; t < c_dif.GetDim4(); ++t)
    for (int k = pmb->ks; k <= pmb->ke; ++k)
      for (int j = pmb->js; j <= pmb->je; ++j) {
        pco->CellVolume(k, j, pmb->is, pmb->ie, vol_);
        for (int i = pmb->is; i <= pmb->ie; ++i)
          if (c_dif(t,k,j,i) > 0.) {
            MaterialPoint p;
            p.ct = t;
            p.x1 = pco->x1f(i) + (1.*rand()/RAND_MAX)*pco->dx1f(i);
            p.x2 = pco->x2f(i) + (1.*rand()/RAND_MAX)*pco->dx2f(i);
            p.x3 = pco->x3f(i) + (1.*rand()/RAND_MAX)*pco->dx3f(i);
            p.v1 = 0.;
            p.v2 = 0.;
            p.v3 = 0.;
            p.mass = c_dif(t,k,j,i)*vol_(i);
            mp.push_back(p);
          }
      }
  c_dif.ZeroClear();
}

void Particles::TimeIntegrate(std::vector<MaterialPoint> &mp, Real time, Real dt)
{
  TranslateEuler(mp, dt);
}

void Particles::WeightedAverage(std::vector<MaterialPoint> &mp_out,
  std::vector<MaterialPoint> const& mp_in, Real ave_wghts[])
{
  size_t psize = mp_out.size();
  for (size_t i = 0; i < psize; ++i) {
    mp_out[i].x1 = ave_wghts[0]*mp_out[i].x1 + ave_wghts[1]*mp_in[i].x1;
    mp_out[i].x2 = ave_wghts[0]*mp_out[i].x2 + ave_wghts[1]*mp_in[i].x2;
    mp_out[i].x3 = ave_wghts[0]*mp_out[i].x3 + ave_wghts[1]*mp_in[i].x3;

    mp_out[i].v1 = ave_wghts[0]*mp_out[i].v1 + ave_wghts[1]*mp_in[i].v1;
    mp_out[i].v2 = ave_wghts[0]*mp_out[i].v2 + ave_wghts[1]*mp_in[i].v2;
    mp_out[i].v3 = ave_wghts[0]*mp_out[i].v3 + ave_wghts[1]*mp_in[i].v3;
  }
}
