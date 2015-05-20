#include "Mesh.h"
#include "Particle.h"
#include "Flow.h"

#ifndef PARTICLECONTAINER_H
#define PARTICLECONTAINER_H

void BuildParticleContainerPlanar( Mesh* Fmesh, long int UmeshN[3], long int UmeshD[3], Particle* &particle, long int numP );
void BuildParticleContainerFullX2( Mesh* mesh, Particle* &pl, long int &numP );
void BuildParticleContainerFull( Mesh* mesh, Particle* &pl, long int &numP );
void BuildParticleContainerOne( Particle* &pl );
void BuildParticleContainerTwo( Particle* &pl );
void BuildParticleContainerCube( Particle* &pl, double bb[], int np );

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

