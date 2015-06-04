#include <Mesh.h>
#include <Particle.h>
#include <Flow.h>

#ifndef PARTICLECONTAINER_H
#define PARTICLECONTAINER_H

#ifndef DO_MPI
#define DO_MPI 0
#endif

void BuildParticleContainerPlanar( Mesh* Fmesh, long int UmeshN[3], double UmeshD[3], Particle* &particle, long int numP, double stepsize );
void BuildParticleContainerFullX2( Mesh* mesh, Particle* &pl, long int &numP, double stepsize );

#if DO_MPI
void BuildParticleContainerFull( Mesh* mesh, Particle* &pl, long int &numP, double stepsize, int rank, int numProcs );
#else
void BuildParticleContainerFull( Mesh* mesh, Particle* &pl, long int &numP, double stepsize );
#endif
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

