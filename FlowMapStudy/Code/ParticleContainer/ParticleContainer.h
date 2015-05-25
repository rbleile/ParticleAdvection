#include <Mesh.h>
#include <Particle.h>
#include <Flow.h>

#ifndef PARTICLECONTAINER_H
#define PARTICLECONTAINER_H

void BuildParticleContainerPlanar( Mesh* Fmesh, long int UmeshN[3], long int UmeshD[3], Particle* &particle, long int numP, double stepsize );
void BuildParticleContainerFullX2( Mesh* mesh, Particle* &pl, long int &numP, double stepsize );
void BuildParticleContainerFull( Mesh* mesh, Particle* &pl, long int &numP, double stepsize );
void BuildParticleContainerOne( Particle* &pl, double stepsize );
void BuildParticleContainerTwo( Particle* &pl, double stepsize );
void BuildParticleContainerCube( Particle* &pl, double bb[], int np, double stepsize );

class ParticleContainer
{
  private:
	long int numParticles;
  public:
	Particle *particle;

	ParticleContainer( Particle *p, long int nump );
	ParticleContainer( Flow *fl, long int nump, double stepsize );

	long int getNumParticles()
	{
		return numParticles;
	}

};

#endif

