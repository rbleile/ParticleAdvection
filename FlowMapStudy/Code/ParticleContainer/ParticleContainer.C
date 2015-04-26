#include "ParticleContainer.h"
#include "Particle.h"
#include "Mesh.h"

ParticleContainer::ParticleContainer( Particle *p, long int nump )
{
	numParticles = nump;
	particle = p;
}

ParticleContainer::ParticleContainer( Flow *fl, long int nump )
{
	numParticles = nump;
	particle = new Particle [numParticles];
	for( long int np=0; np < numParticles; np++ )
	{
		particle[np].setPoint( fl[np].out.x, fl[np].out.y, fl[np].out.z );
		particle[np].setStepSize( 0.001 );
	}
}

void BuildParticleContainerFullX2( Mesh* mesh, Particle *pl, long int &numP )
{
	long int nx = mesh->nx*2;
	long int ny = mesh->ny*2;
	long int nz = mesh->nz*2;

	double x0 = mesh->x0;
	double y0 = mesh->y0;
	double z0 = mesh->z0;

	double dx = mesh->dx / 2.0;
	double dy = mesh->dy / 2.0;
	double dz = mesh->dz / 2.0;

	numP = nx * ny * nz;	

	pl = new Particle [ numP ];

	for( int z = 0; z < nz; z++ )
	{
		for( int y = 0; y < ny; y++ )
		{
			for( int x = 0; x < nx; x++ )
			{
				long int id = x + nx * ( y + ( ny * z ) );

				pl[id].x = x0 + x*dx;
				pl[id].y = y0 + y*dy;
				pl[id].z = z0 + z*dz;
				pl[id].setStepSize( STEPSIZE );

			}
		}
	}

}


void BuildParticleContainerFull( Mesh* mesh, Particle *pl, long int &numP )
{
	long int nx = mesh->nx;
	long int ny = mesh->ny;
	long int nz = mesh->nz;

	double x0 = mesh->x0;
	double y0 = mesh->y0;
	double z0 = mesh->z0;

	double dx = mesh->dx;
	double dy = mesh->dy;
	double dz = mesh->dz;

	numP = nx * ny * nz;	

	pl = new Particle [ numP ];

	for( int z = 0; z < nz; z++ )
	{
		for( int y = 0; y < ny; y++ )
		{
			for( int x = 0; x < nx; x++ )
			{
				long int id = x + nx * ( y + ( ny * z ) );

				pl[id].x = x0 + x*dx;
				pl[id].y = y0 + y*dy;
				pl[id].z = z0 + z*dz;
				pl[id].setStepSize( STEPSIZE );

			}
		}
	}

}

void BuildParticleContainerOne( Particle* &pl )
{
	pl = new Particle [1];
	pl[0].x = 1.0;
	pl[0].y = 1.0;
	pl[0].z = 1.0;
	pl[0].setStepSize( STEPSIZE );
}

void BuildParticleContainerTwo( Particle* &pl )
{
	pl = new Particle [2];
	pl[0].x = 1.0;
	pl[0].y = 1.0;
	pl[0].z = 1.0;
	pl[0].setStepSize( STEPSIZE );

	pl[1].x = -1.0;
	pl[1].y = -1.0;
	pl[1].z = -1.0;
	pl[1].setStepSize( STEPSIZE );
}

