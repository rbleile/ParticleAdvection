#include <cstdio>
#include <iostream>

#include "AdvectParticles.h"

#include <sys/time.h>
#define GET_TIME(now){ \
    struct timeval t; \
    gettimeofday(&t, NULL); \
    now = t.tv_sec + t.tv_usec/1000000.0; \
}

using std::cerr;
using std::cout;
using std::endl;

long int GetParticleCellID( Mesh* FineMesh, DomainMesh *UberMesh, Particle* part )
{
	double vel[3];
	FineMesh->getVelocity( *part, vel );

	double xmax = FineMesh->x0 + FineMesh->dx*(FineMesh->nx-1);
	double ymax = FineMesh->y0 + FineMesh->dy*(FineMesh->ny-1);
	double zmax = FineMesh->z0 + FineMesh->dz*(FineMesh->nz-1);

	Point np( part->x + vel[0]*STEPSIZE, part->y + vel[1]*STEPSIZE, part->z + vel[2]*STEPSIZE );

	long int cellID = -1;

	if( np.x > xmax || np.x < FineMesh->x0 
	 || np.y > ymax || np.y < FineMesh->y0
	 || np.z > zmax || np.z < FineMesh->z0 )
	{
		return cellID;
	}
	else
	{
		cellID = UberMesh->getCellID( np.x, np.y, np.z );

		if( np.x == xmax  )
		{
			cellID -= 1;
		}
		if( np.y == ymax  )
		{
			cellID -= UberMesh->nx;
		}
		if( np.z == zmax  )
		{
			cellID -= UberMesh->nx*UberMesh->ny;
		}		
	}
	return cellID;
}

void ComputeFaceID( long int cellID, DomainMesh *UberMesh, Particle* part )
{
	double bb[6];
	UberMesh->getCellBounds( cellID, bb );
	part->setFID( bb );
}

bool AcceptableFlow( long int cellID, DomainMesh *UberMesh, Mesh* FineMesh, Particle* part )
{

	int FID = part->FID;

	if( FID < 0 )
	{
		return 0;
	}
	
	long int cell_id = -1;
	int d;

	if( FID < 6 )
	{
		d = (( FID == 1 || FID == 3) ? 0 : ( ( FID == 0 || FID == 2  ) ? 1  : 2 ) );

		long int Rcdim[3];
		if     ( d == 0 ) { Rcdim[0] = UberMesh->nx-1;	Rcdim[1] = FineMesh->ny-1;	Rcdim[2] = FineMesh->nz-1; }
		else if( d == 1 ) { Rcdim[0] = FineMesh->nx-1;	Rcdim[1] = UberMesh->ny-1;	Rcdim[2] = FineMesh->nz-1; }
		else if( d == 2 ) { Rcdim[0] = FineMesh->nx-1;	Rcdim[1] = FineMesh->ny-1;	Rcdim[2] = UberMesh->nz-1;   }

		long int Uids[3];
		long int Fids[3];

		UberMesh->D1to3P( cellID, Uids );
		FineMesh->getLogicalCellID( part->x, part->y, part->z, Fids );

		long int x_id, y_id, z_id;

		if      ( d == 0 ){ x_id = Uids[0]; y_id = Fids[1]; z_id = Fids[2]; }
		else if ( d == 1 ){ x_id = Fids[0]; y_id = Uids[1]; z_id = Fids[2]; }
		else if ( d == 2 ){ x_id = Fids[0]; y_id = Fids[1]; z_id = Uids[2]; }

		cell_id = x_id + Rcdim[0] * ( y_id + ( Rcdim[1] * z_id ) );
	}
	else if ( FID >= 6 )
	{
		cerr << "FID: EDGE CASE HIT: " << FID << "\t This case is not yet handled" << endl;
		return 0;
	}	

	return UberMesh->AcceptableFlow[d][ cell_id ];

}

int AdvectParticleOnFlow( int cellID, Mesh *FineMesh, DomainMesh *UberMesh, Particle* part )
{
	int FID = part->FID;

	if( FID < 0 )
	{
		return 2;
	}
	
	long int base_id;
	int d;

	if( FID < 6 )
	{
		d = (( FID == 1 || FID == 3) ? 0 : ( ( FID == 0 || FID == 2  ) ? 1  : 2 ) );

		long int Rdim[3];
		if     ( d == 0 ) { Rdim[0] = UberMesh->nx;	Rdim[1] = FineMesh->ny;	Rdim[2] = FineMesh->nz; }
		else if( d == 1 ) { Rdim[0] = FineMesh->nx;	Rdim[1] = UberMesh->ny;	Rdim[2] = FineMesh->nz; }
		else if( d == 2 ) { Rdim[0] = FineMesh->nx;	Rdim[1] = FineMesh->ny;	Rdim[2] = UberMesh->nz; }

		long int xstep = 1;
		long int ystep = Rdim[0];
		long int zstep = Rdim[0]*Rdim[1];

		long int Rneighbors[3];
		if( d == 0 )     { Rneighbors[0] = ystep;	Rneighbors[1] = zstep;	Rneighbors[2] = ystep+zstep; }
		else if( d == 1 ){ Rneighbors[0] = zstep;	Rneighbors[1] = xstep;	Rneighbors[2] = zstep+xstep; }
		else if( d == 2 ){ Rneighbors[0] = xstep;	Rneighbors[1] = ystep;	Rneighbors[2] = xstep+ystep; }

		long int Uids[3];
		long int Fids[3];

		UberMesh->D1to3P( cellID, Uids );
		FineMesh->getLogicalCellID( part->x, part->y, part->z, Fids );

		long int x_id, y_id, z_id;

		if      ( d == 0 ){ x_id = Uids[0]; y_id = Fids[1]; z_id = Fids[2]; }
		else if ( d == 1 ){ x_id = Fids[0]; y_id = Uids[1]; z_id = Fids[2]; }
		else if ( d == 2 ){ x_id = Fids[0]; y_id = Fids[1]; z_id = Uids[2]; }

		base_id = x_id + Rdim[0] * ( y_id + ( Rdim[1] * z_id ) );

		long int flowIDs[4] = { base_id, base_id + Rneighbors[0], base_id + Rneighbors[1], base_id + Rneighbors[2] };

		double x = part->x; double y = part->y; double z = part->z; double t = part->t;
		double fracD1, fracD2;

		Flow* flows = UberMesh->flowField[d];

		Point p0 = flows[ flowIDs[0] ].in;
		Point p1 = flows[ flowIDs[1] ].in;
		Point p2 = flows[ flowIDs[2] ].in;
		Point p3 = flows[ flowIDs[3] ].in;

		if      ( d == 0 ){ fracD1 = ( y - p0.y ) / ( p1.y - p0.y ); fracD2 = ( z - p0.z )/( p2.z - p0.z ); }
		else if ( d == 1 ){ fracD1 = ( z - p0.z ) / ( p1.z - p0.z ); fracD2 = ( x - p0.x )/( p2.x - p0.x ); }
		else if ( d == 2 ){ fracD1 = ( x - p0.x ) / ( p1.x - p0.x ); fracD2 = ( y - p0.y )/( p2.y - p0.y ); }

		p0 = flows[ flowIDs[0] ].out;
		p1 = flows[ flowIDs[1] ].out;
		p2 = flows[ flowIDs[2] ].out;
		p3 = flows[ flowIDs[3] ].out;
	
		double X01 = p0.x + fracD1 * ( p1.x - p0.x ); 
		double X23 = p0.x + fracD1 * ( p3.x - p2.x );

		double Y01 = p0.y + fracD1 * ( p1.y - p0.y ); 
		double Y23 = p2.y + fracD1 * ( p3.y - p2.y );

		double Z01 = p0.z + fracD1 * ( p1.z - p0.z ); 
		double Z23 = p2.z + fracD1 * ( p3.z - p2.z );

		double T01 = p0.t + fracD1 * ( p1.t - p0.t ); 
		double T23 = p2.t + fracD1 * ( p3.t - p2.t );

		double xx = x;
		double yy = y;
		double zz = z;
		double tt = t;
		
		x = X01 + fracD2 * ( X23-X01 );
		y = Y01 + fracD2 * ( Y23-Y01 );
		z = Z01 + fracD2 * ( Z23-Z01 );
		t = T01 + fracD2 * ( T23-T01 );
	
		if( tt == t )
		{
			std::cerr << "dt = 0 " << std::endl;
			return 0;
		}

		part->x = x;
		part->y = y;
		part->z = z;
		part->t += t;

		return 1;

	}
	else if ( FID >= 6 )
	{
		cerr << "FID: EDGE CASE HIT: " << FID << "\t This case is not yet handled" << endl;
		return 0;
	}	

}

void AdvectParticleList( Mesh *FineMesh, DomainMesh *UberMesh, ParticleContainer* advectList, double endtime, double* MBB )
{

	int numParticles = advectList->getNumParticles();

	for( long int i = 0; i < numParticles; i++ ){

		int status = 1;
		Particle &particle = advectList->particle[i];

		long int cellID = GetParticleCellID( FineMesh, UberMesh, &particle );
		double cbb[6];

		if( cellID == -1 ) continue; //Particle Seeded outside of Mesh ( Skip it )

		UberMesh->getCellBounds( cellID, cbb );

		int onCellFace = FineMesh->onBoundary( particle, cbb );	

		if( onCellFace ) ComputeFaceID( cellID, UberMesh, &particle );

		int count = 0;
		double start_time, end_time;
		double AdvectionTime = 0;
		double EulerTime = 0;

		while( status )
		{	
			if( onCellFace )
			{
				ComputeFaceID( cellID, UberMesh, &particle );
				int canFlow = AcceptableFlow( cellID, UberMesh, FineMesh, &particle );
				if( canFlow )
				{
					GET_TIME( start_time );
					AdvectParticleOnFlow( cellID, FineMesh, UberMesh, &particle );
					GET_TIME( end_time );
					AdvectionTime += end_time-start_time;
					count++;
				}
				else
				{
					double t = particle.t;
					GET_TIME( start_time );
					status = UberMesh->EulerCellAdvection( cellID, endtime, MBB, particle ); 
					GET_TIME( end_time );
					EulerTime += end_time-start_time;
					double tt = particle.t;
					count += (int)(tt-t)/STEPSIZE;
				}
			}
			else
			{
				double t = particle.t;
				GET_TIME( start_time );
				status = UberMesh->EulerCellAdvection( cellID, endtime, MBB, particle ); 
				GET_TIME( end_time );
				EulerTime += end_time-start_time;
				double tt = particle.t;
				count += (int)(tt-t)/STEPSIZE;
			}

			if( status == 2 )
			{
				cellID = GetParticleCellID( FineMesh, UberMesh, &particle ); 
				ComputeFaceID( cellID, UberMesh, &particle );
				onCellFace = 1;
			}


			status = ( particle.t >= endtime ) ? 0 : status;
		}	

		cout << "Advection Time: " << AdvectionTime << endl;
		cout << "Euler Time: " << EulerTime << endl;
		cout << "EndTime: " << particle.t << endl;
		cout << "Count:       " << count << endl;
		cout << "Total Euler: " << endtime/STEPSIZE << endl;

		cerr << "Finished Particle at:" << endl;
		cerr << particle.x << " " << particle.y << " " << particle.z  << " " << particle.t << endl;
		
	}
}
