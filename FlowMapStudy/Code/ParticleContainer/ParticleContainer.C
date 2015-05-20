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

void BuildParticleContainerPlanar( Mesh* Fmesh, long int UmeshN[3], long int UmeshD[3], Particle* &particle, long int numP )
{

	particle = new Particle [ numP ];

	double Rpos[3] = { Fmesh->x0, Fmesh->y0, Fmesh->z0 };
	for( int d = 0; d < 3; d++ )
	{
		long int Rdims[3];
		if( d == 0 )     { Rdims[0] = UmeshN[0];	Rdims[1] = Fmesh->ny;	Rdims[2] = Fmesh->nz; }
		else if( d == 1 ){ Rdims[0] = Fmesh->nx;	Rdims[1] = UmeshN[1];	Rdims[2] = Fmesh->nz; }
		else if( d == 2 ){ Rdims[0] = Fmesh->nx;	Rdims[1] = Fmesh->ny;	Rdims[2] = UmeshN[2]; }


		double Rdel[3];
		if( d == 0 )     { Rdel[0] = UmeshD[0];	Rdel[1] = Fmesh->dy;	Rdel[2] = Fmesh->dz; }
		else if( d == 1 ){ Rdel[0] = Fmesh->dx;	Rdel[1] = UmeshD[1];	Rdel[2] = Fmesh->dz; }
		else if( d == 2 ){ Rdel[0] = Fmesh->dx;	Rdel[1] = Fmesh->dy;	Rdel[2] = UmeshD[2]; }

		for( long int z_id = 0; z_id < Rdims[2]; z_id++ )
		{	
			for( long int y_id = 0; y_id < Rdims[1]; y_id++ )
			{	
				for( long int x_id = 0; x_id < Rdims[0]; x_id++ )
				{	
	
					long int id = x_id + Rdims[0] * ( y_id + Rdims[1] * z_id ) + Rdims[0]*Rdims[1]*Rdims[2]*d;

					double x = Rpos[0] + Rdel[0]*x_id;
					double y = Rpos[1] + Rdel[1]*y_id;
					double z = Rpos[2] + Rdel[2]*z_id;

					particle[id].x = x;
					particle[id].y = y;
					particle[id].z = z;
					particle[id].t = 0.0;
					particle[id].setStepSize( STEPSIZE );
				}
			}
		}
	}
}

void BuildParticleContainerFullX2( Mesh* mesh, Particle* &pl, long int &numP )
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


void BuildParticleContainerFull( Mesh* mesh, Particle* &pl, long int &numP )
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

	cerr << "IN: " << pl[0].x << " " << pl[0].y << " " << pl[0].z << endl;

}

void BuildParticleContainerOne( Particle* &pl )
{
	pl = new Particle [1];
	pl[0].x = 0.646207;
	pl[0].y = 1.55423;
	pl[0].z = -0.0348686;
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

void BuildParticleContainerCube( Particle* &pl, double bb[], int np )
{
	pl = new Particle [np*np*np];

	double dx = (bb[1]-bb[0]) / np;
	double dy = (bb[3]-bb[2]) / np;
	double dz = (bb[5]-bb[4]) / np;	

	for( int z = 0; z < np; z++ ){
	for( int y = 0; y < np; y++ ){
	for( int x = 0; x < np; x++ )
	{
		int i = x + y*np + z*np*np;
		pl[i].x = bb[0] + dx*x;
		pl[i].y = bb[2] + dy*y;
		pl[i].z = bb[4] + dz*z;
		pl[i].setStepSize( STEPSIZE );
	}}}
}

