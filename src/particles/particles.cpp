/** @file particles.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Sunday May 30, 2021 13:10:12 PDT
 * @bug No known bugs.
 */

// C++ headers
#include <ctime>
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "particles.hpp"
#include "particle_buffer.hpp"

Particles::Particles(MeshBlock *pmb, ParameterInput *pin):
  pmy_block(pmb), myname("HEAD"), prev(nullptr), next(nullptr),
  seeds_per_cell_(1), nmax_per_cell_(1), density_floor_(0), has_gravity_(false)
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
    } else if (std::strncmp(p, "scp", 3) == 0) {
      AddParticles(SimpleCloudParticles(pmb, pin, name));
    } else {
      msg << "### FATAL ERROR in function Particles::Particles"
          << std::endl << "Particles '" << p << "' "
          << "unrecognized." << std::endl;
      ATHENA_ERROR(msg);
    }
    p = std::strtok(NULL, " ,");
  }
}

// constructor, initializes data structure and parameters
Particles::Particles(MeshBlock *pmb, ParameterInput *pin, std::string name, int nct):
  pmy_block(pmb), myname(name), prev(nullptr), next(nullptr)
{
  ppb = new ParticleBuffer(this);
  int nc1 = pmb->ncells1, nc2 = pmb->ncells2, nc3 = pmb->ncells3;

  xface_.resize(nc3+nc2+nc1+3);
  for (int k = 0; k <= nc3; ++k)
    xface_[k] = pmb->pcoord->x3f(k);
  for (int j = 0; j <= nc2; ++j)
    xface_[nc3+1+j] = pmb->pcoord->x2f(j);
  for (int i = 0; i <= nc1; ++i)
    xface_[nc3+nc2+2+i] = pmb->pcoord->x1f(i);

  xcenter_.resize(nc3+nc2+nc1);
  for (int k = 0; k < nc3; ++k)
    xcenter_[k] = pmb->pcoord->x3v(k);
  for (int j = 0; j < nc2; ++j)
    xcenter_[nc3+j] = pmb->pcoord->x2v(j);
  for (int i = 0; i < nc1; ++i)
    xcenter_[nc3+nc2+i] = pmb->pcoord->x1v(i);

  dims_.resize(3);
  dims_[0] = nc3;
  dims_[1] = nc2;
  dims_[2] = nc1;

  c.NewAthenaArray(nct, nc3, nc2, nc1);
  c.ZeroClear();
  c1_.NewAthenaArray(nct, nc3, nc2, nc1);
  c1_.ZeroClear();
  pcell_.NewAthenaArray(nct, nc3, nc2, nc1);

  seeds_per_cell_ = pin->GetOrAddInteger("particles", name + ".seeds_per_cell", 5);
  nmax_per_cell_ = pin->GetOrAddInteger("particles", name + ".nmax_per_cell",
    5*seeds_per_cell_);
  density_floor_ = pin->GetOrAddReal("particles", name + ".dfloor", 1.E-10);
  has_gravity_ = false;

  std::stringstream msg;
  if (nmax_per_cell_ < seeds_per_cell_) {
    msg << "### FATAL ERROR in Particles::Particles"
        << "Maximum particles per cell: " << nmax_per_cell_
        << " is less than "
        << "seed particles per cell: " << seeds_per_cell_ << std::endl;
    ATHENA_ERROR(msg);
  }
}

// destructor
Particles::~Particles()
{
  if (prev != nullptr) prev->next = next;
  if (next != nullptr) next->prev = prev;
  delete ppb;
}

Particles::Particles(Particles const& other):
  c(other.c), c1_(other.c1_), pcell_(other.pcell_)
{
  if (this == &other) return;
  *this = other;
  ppb = new ParticleBuffer(this);
}

/*Particles& Particles::operator=(Particles const& other)
{
  pmy_block = other.pmy_block;
  myname = other.myname;
  prev = other.prev;
  next = other.next;
  // assignment operator of AthenaArray does not allocate memory
  c = other.c;
  mp = other.mp;
  mp1 = other.mp1;

  xface_ = other.xface_;
  xcenter_ = other.xcenter_;
  dims_ = other.dims_;
  cnames_ = other.cnames_;
  available_ids_ = other.available_ids_;
  cc_ = other.cc_;
  mu_ = other.mu_;
  c1_ = other.c1_;
  pcell_ = other.pcell_;
  seeds_per_cell_ = other.seeds_per_cell_;
  nmax_per_cell_ = other.nmax_per_cell_;
  density_floor_ = other.density_floor_;
  has_gravity_ = other.has_gravity_;

  ppb = new ParticleBuffer(this);

  return *this;
}*/

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

void Particles::Initialize()
{
  Particles *p = this;
  while (p != nullptr) {
    Particulate(p->mp, p->c);
    p = p->next;
  }
}

void Particles::TimeIntegrate(std::vector<MaterialPoint> &mp, Real time, Real dt)
{
  for (std::vector<MaterialPoint>::iterator it = mp.begin(); it != mp.end(); ++it) {
    it->x1 += it->v1*dt;
    it->x2 += it->v2*dt;
    it->x3 += it->v3*dt;
  }
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
