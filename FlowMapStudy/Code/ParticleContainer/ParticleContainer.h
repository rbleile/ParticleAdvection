#include "Mesh.h"
#include "Particle.h"
#include "Flow.h"

#ifndef PARTICLECONTAINER_H
#define PARTICLECONTAINER_H

void BuildParticleContainerFullX2( Mesh* mesh, Particle *pl, long int &numP );
void BuildParticleContainerFull( Mesh* mesh, Particle *pl, long int &numP );
void BuildParticleContainerOne( Particle* &pl );

class ParticleContainer
{
  private:
	long int numParticles;
  public:
	Particle *particle;

	ParticleContainer( Particle *p, long int nump );
	ParticleContainer( Flow *fl, long int nump );

	long int getNumParticles()
	{
		return numParticles;
	}

};

#endif

