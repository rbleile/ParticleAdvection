#include "MeshParticle.h"

#ifndef MESH_PARTICLE_LIST
#define MESH_PARTICLE_LIST

class ParticleList
{
  private:
	long int numParticles;
  public:
	Particle *particle;

	ParticleList( Particle *p, long int nump );
	ParticleList( FMFlow *fl, long int nump );

	long int getNumParticles()
	{
		return numParticles;
	}

};

#endif

