#include "MeshParticleList.h"

ParticleList::ParticleList( Particle *p, long int nump )
{
	numParticles = nump;
	particle = p;
}

ParticleList::ParticleList( FMFlow *fl, long int nump )
{
	numParticles = nump;
	particle = new Particle [numParticles];
	for( long int np=0; np < numParticles; np++ )
	{
		particle[np].setPoint( fl[np].out.x, fl[np].out.y, fl[np].out.z );
		particle[np].setStepSize( 0.001 );
	}
}

