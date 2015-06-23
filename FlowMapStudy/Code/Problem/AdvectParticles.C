#include <cstdio>
#include <iostream>
#include <cmath>
#include <sstream>
#include <cstring>

#include <AdvectParticles.h>

#include <sys/time.h>
#define GET_TIME(now){ \
    struct timeval t; \
    gettimeofday(&t, NULL); \
    now = t.tv_sec + t.tv_usec/1000000.0; \
}

#ifndef doCopy
#define doCopy 0
#endif

#ifndef DO_MPI
#define DO_MPI 0
#endif

using std::cerr;
using std::cout;
using std::endl;
using std::sqrt;
using std::pow;
using std::fabs;
using std::min;
using std::max;
using std::string;
using std::stringstream;

int RANK = 0;
int ITERPRIME = 0;

const int FacesToFaces[6]= { 2, 3, 0, 1, 5, 4 };
inline int checkBoundary( Particle particle, double* bbox )
{

    if( particle.x > bbox[0] && particle.x < bbox[1] &&
        particle.y > bbox[2] && particle.y < bbox[3] &&
        particle.z > bbox[4] && particle.z < bbox[5] )
    {
        return 1;
    }
    else
    {
        return 0;
    }
}
long int GetParticleCellID( Mesh* FineMesh, DomainMesh *UberMesh, Particle* part, double endTime, double* MBB)
{
	
	if( !checkBoundary( *part, MBB ))
	{
		return -1;
	}
	
	double xmax = FineMesh->x0 + FineMesh->dx*(FineMesh->nx-1);
	double ymax = FineMesh->y0 + FineMesh->dy*(FineMesh->ny-1);
	double zmax = FineMesh->z0 + FineMesh->dz*(FineMesh->nz-1);

	Particle np( part );

	//if(RANK == 0 && ITERPRIME == 13631 ) cerr << "CELLID FUNC 4: " << np.x << " " << np.y << " " << np.z << " " << np.t << endl;
	int stat = FineMesh->RK4( MBB, MBB, endTime, np );
	//if(RANK == 0 && ITERPRIME == 13631 ) cerr << "stat: " << stat << endl;
	if( stat == 0 ){ return -1; }

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

		UberMesh->D1to3C( cellID, Uids );
		FineMesh->getLogicalCellID( part->x, part->y, part->z, Fids );

		Uids[0] = ( FID == 1 ) ? Uids[0]+1 : Uids[0];
		Uids[1] = ( FID == 2 ) ? Uids[1]+1 : Uids[1];
		Uids[2] = ( FID == 5 ) ? Uids[2]+1 : Uids[2];

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

int AdvectParticleOnFlow( long int &cellID, Mesh *FineMesh, DomainMesh *UberMesh, Particle* part, double endTime, double* MBB )
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


#if doCopy

					fprintf( stderr, "doCopy\n", i );
					Particle copy;
					copy.x = part->x;
					copy.y = part->y;
					copy.z = part->z;
					copy.t = part->t;
					copy.setStepSize( part->getStepSize() );
					UberMesh->EulerCellAdvection( cellID, endTime, MBB, copy ); 

#endif			 


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

		UberMesh->D1to3C( cellID, Uids );
		FineMesh->getLogicalCellID( part->x, part->y, part->z, Fids );

		Uids[0] = ( FID == 1 ) ? Uids[0]+1 : Uids[0];
		Uids[1] = ( FID == 2 ) ? Uids[1]+1 : Uids[1];
		Uids[2] = ( FID == 5 ) ? Uids[2]+1 : Uids[2];

		long int x_id, y_id, z_id;

		if      ( d == 0 ){ x_id = Uids[0]; y_id = Fids[1]; z_id = Fids[2]; }
		else if ( d == 1 ){ x_id = Fids[0]; y_id = Uids[1]; z_id = Fids[2]; }
		else if ( d == 2 ){ x_id = Fids[0]; y_id = Fids[1]; z_id = Uids[2]; }

		base_id = x_id + Rdim[0] * ( y_id + ( Rdim[1] * z_id ) );
		
		long int flowIDs[4] = { base_id, base_id + Rneighbors[0], base_id + Rneighbors[1], base_id + Rneighbors[2] };

		double x = part->x; double y = part->y; double z = part->z; double t = part->t;
		double fracD1, fracD2;

		Flow* flows = UberMesh->flowField[d];

		int ordering[4];
		Point pi[4]= {	flows[ flowIDs[0] ].in,
						flows[ flowIDs[1] ].in, 
						flows[ flowIDs[2] ].in, 
						flows[ flowIDs[3] ].in  };
	
		double minX = min( pi[0].x, min( pi[1].x, min( pi[2].x, pi[3].x ) ) );
		double minY = min( pi[0].y, min( pi[1].y, min( pi[2].y, pi[3].y ) ) );
		double minZ = min( pi[0].z, min( pi[1].z, min( pi[2].z, pi[3].z ) ) );

		double maxX = max( pi[0].x, max( pi[1].x, max( pi[2].x, pi[3].x ) ) );
		double maxY = max( pi[0].y, max( pi[1].y, max( pi[2].y, pi[3].y ) ) );
		double maxZ = max( pi[0].z, max( pi[1].z, max( pi[2].z, pi[3].z ) ) );

		double f11,f12,f21,f22;

		if( d == 0 )
		{
			//Y Z ordering	
			double denom = (maxY-minY)*(maxZ-minZ);

			f11 = (( maxY - y ) * ( maxZ - z )) / ( denom );
			f21 = (( y - minY ) * ( maxZ - z )) / ( denom );
			f12 = (( maxY - y ) * ( z - minZ )) / ( denom );
			f22 = (( y - minY ) * ( z - minZ )) / ( denom );

			for( int i = 0; i < 4; i++ )
			{
				if     ( pi[i].y == minY && pi[i].z == minZ )
				{
					ordering[0] = i;
				}
				else if( pi[i].y == maxY && pi[i].z == minZ )
				{
					ordering[1] = i;
				}
				else if( pi[i].y == minY && pi[i].z == maxZ )
				{
					ordering[2] = i;
				}
				else if( pi[i].y == maxY && pi[i].z == maxZ )
				{
					ordering[3] = i;
				}
			}

		}
		else if( d == 1 )
		{
		//Z X ordering
			double denom = (maxZ-minZ)*(maxX-minX);

			f11 = (( maxZ - z ) * ( maxX - x )) / ( denom );
			f21 = (( z - minZ ) * ( maxX - x )) / ( denom );
			f12 = (( maxZ - z ) * ( x - minX )) / ( denom );
			f22 = (( z - minZ ) * ( x - minX )) / ( denom );

			for( int i = 0; i < 4; i++ )
			{
				if     ( pi[i].z == minZ && pi[i].x == minX )
				{
					ordering[0] = i;
				}
				else if( pi[i].z == maxZ && pi[i].x == minX )
				{
					ordering[1] = i;
				}
				else if( pi[i].z == minZ && pi[i].x == maxX )
				{
					ordering[2] = i;
				}
				else if( pi[i].z == maxZ && pi[i].x == maxX )
				{
					ordering[3] = i;
				}
			}
		}
		else if( d == 2 )
		{
		//X Y ordering
			double denom = (maxX-minX)*(maxY-minY);

			f11 = (( maxX - x ) * ( maxY - y )) / ( denom );
			f21 = (( x - minX ) * ( maxY - y )) / ( denom );
			f12 = (( maxX - x ) * ( y - minY )) / ( denom );
			f22 = (( x - minX ) * ( y - minY )) / ( denom );

			for( int i = 0; i < 4; i++ )
			{
				if     ( pi[i].x == minX && pi[i].y == minY )
				{
					ordering[0] = i;
				}
				else if( pi[i].x == maxX && pi[i].y == minY )
				{
					ordering[1] = i;
				}
				else if( pi[i].x == minX && pi[i].y == maxY )
				{
					ordering[2] = i;
				}
				else if( pi[i].x == maxX && pi[i].y == maxY )
				{
					ordering[3] = i;
				}
			}
		}

			
//same for all
			Point po[4]= {	flows[ flowIDs[0] ].out,
							flows[ flowIDs[1] ].out, 
							flows[ flowIDs[2] ].out, 
							flows[ flowIDs[3] ].out  };
			double tt = t;

			x = f11*po[ordering[0]].x + f21*po[ordering[1]].x + f12*po[ordering[2]].x + f22*po[ordering[3]].x;
			y = f11*po[ordering[0]].y + f21*po[ordering[1]].y + f12*po[ordering[2]].y + f22*po[ordering[3]].y;
			z = f11*po[ordering[0]].z + f21*po[ordering[1]].z + f12*po[ordering[2]].z + f22*po[ordering[3]].z;
			t = f11*po[ordering[0]].t + f21*po[ordering[1]].t + f12*po[ordering[2]].t + f22*po[ordering[3]].t;
	
		//If no time change in advection
		if( tt == t )
		{
			std::cerr << "dt = 0 " << std::endl;
			return 0;
		}

		//If the interpolation will put us past the end time and half of the time difference is greater then end time just euler to finish
		if( (tt + t > endTime) && (tt + (t/2.0) > endTime) )
		{
			UberMesh->EulerCellAdvection( cellID, endTime, MBB, *part ); 
			return 0;
		}

		part->x = x;
		part->y = y;
		part->z = z;
		part->t += t;

		// If the interpolation step has put us past the end time then euler backwards to get to end time.
		if( part->t > endTime )
		{	
			t = part->t;
			
			UberMesh->ReverseEulerCellAdvection( cellID, endTime, MBB, *part );
			return 0;
		}


		//Compute the face we  landed on
		ComputeFaceID( cellID, UberMesh, part );

#if doCopy

		long int dcell_id = x_id + (Rdim[0]-1) * ( y_id + ( (Rdim[1]-1) * z_id ) );
		double distance = fabs( sqrt( pow( (part->x-copy.x), 2 ) + pow( (part->y-copy.y), 2 ) + pow( (part->z-copy.z),2) ) );

		double* max_diff = UberMesh->max_diff[d];

		#pragma omp critical
		{
			max_diff[dcell_id] = ( max_diff[dcell_id] < distance  ) ? distance : max_diff[dcell_id];
		}

#endif
	

		return 1;

	}
	else if ( FID >= 6 )
	{
		cerr << "FID: EDGE CASE HIT: " << FID << "\t This case is not yet handled" << endl;
		return 0;
	}	

}

void AdvectParticleList( Mesh* FineMesh, ParticleContainer* advectList, double endtime, double* MBB )
{
	long int numParticles = advectList->getNumParticles();

	int chunkSize = (int)((numParticles/8.0)*0.1);
	if( chunkSize < 1) chunkSize = 1;
	#pragma omp parallel for schedule (dynamic,chunkSize)
	for( long int i = 0; i < numParticles; i++ )
	{
		int count = 0;
		int status = 1;
		Particle &particle = advectList->particle[i];
	
		if( !checkBoundary( particle, MBB ))
		{
			continue;
		}

		while( status == 1 )
		{
			status = FineMesh->RK4( MBB, MBB, endtime, particle ); 
		}
	
	}
}

#if DO_MPI
void AdvectParticleList( Mesh *FineMesh, DomainMesh *UberMesh, ParticleContainer* advectList, double endtime, double* MBB, int &total_lagrange, int &total_euler, int rank, int numProcs)
#else
void AdvectParticleList( Mesh *FineMesh, DomainMesh *UberMesh, ParticleContainer* advectList, double endtime, double* MBB, int &total_lagrange, int &total_euler)
#endif
{

	long int numParticles = advectList->getNumParticles();

	long int nonZeroParticle = 0;

	int totalL = 0;
	int totalE = 0;

	int chunkSize = (int)((numParticles/8.0)*0.1);

	if( chunkSize < 1) chunkSize = 1;
	#pragma omp parallel for shared ( nonZeroParticle, totalL, totalE ) schedule(dynamic,chunkSize)
	for( long int i = 0; i < numParticles; i++ )
	{

		Particle &particle = advectList->particle[i];	

		long int cellID = GetParticleCellID( FineMesh, UberMesh, &particle, endtime, MBB );

		if( cellID == -1 ){
			continue; //Particle Seeded outside of Mesh ( Skip it )
		}

		int status = 1;
		double cbb[6];

		double stepsize = particle.getStepSize();

		UberMesh->getCellBounds( cellID, cbb );

		int onCellFace = FineMesh->onBoundary( particle, cbb );	

		if( onCellFace ) ComputeFaceID( cellID, UberMesh, &particle );

		int count = 0;
		long int total = 0;
		long int advected = 0;
		long int euler = 0;
		double start_time, end_time;


		while( status && count < (int)( endtime/stepsize ) )
		{	
			total++;

			if( onCellFace )
			{
				ComputeFaceID( cellID, UberMesh, &particle );
				int canFlow = AcceptableFlow( cellID, UberMesh, FineMesh, &particle );

				if( canFlow )
				{

					advected++;

					double t = particle.t;
					status = AdvectParticleOnFlow( cellID, FineMesh, UberMesh, &particle, endtime, MBB );
					double tt = particle.t;

					if( tt-t==0.0 )
					{ 
						FineMesh->RK4( MBB, MBB, endtime, particle );	
						tt = particle.t;
					} 
					count++;

				}
				else
				{
					
					double t = particle.t;
					status = UberMesh->EulerCellAdvection( cellID, endtime, MBB, particle ); 
					double tt = particle.t;

					if( tt-t==0.0 )
					{ 
						FineMesh->RK4( MBB, MBB, endtime, particle );	
						tt = particle.t;

					} 

					count += (int)((tt-t)/stepsize);
					euler++;

					
				}
			}
			else
			{
				
				double t = particle.t;
				status = UberMesh->EulerCellAdvection( cellID, endtime, MBB, particle ); 
				double tt = particle.t;
				if( tt-t==0.0 )
				{ 
					FineMesh->RK4( MBB, MBB, endtime, particle );	
					tt = particle.t;

				} 
				count += (int)((tt-t)/stepsize);
				euler++;
					
			}
						
			cellID = GetParticleCellID( FineMesh, UberMesh, &particle, endtime, MBB ); 
			if( cellID == -1 )
			{
				status = 0;
			}
			else{
				ComputeFaceID( cellID, UberMesh, &particle );
				onCellFace = 1;
			}

			status = ( particle.t >= endtime ) ? 0 : status;
		}


		#pragma omp critical
		{
			nonZeroParticle++;
		}	
		#pragma omp critical
		{	
			totalL += advected;
		}
		#pragma omp critical
		{
			totalE += euler;
		}
	}

	total_lagrange += totalL;
	total_euler += totalE;

/*
	cerr << "Total Lagrangian Steps: " << total_lagrange << endl;
	cerr << "Total Eulerian Steps:   " << total_euler << endl;
	cerr << "Average Lagrange PP:    " << (double)total_lagrange / (double)numParticles << endl;
	cerr << "Average Eulerian PP:    " << (double)total_euler / (double)numParticles << endl;
	cerr << "Percent Lagrange Steps: " << ((double)total_lagrange / (double)( total_lagrange + total_euler ) )* 100.0 << endl;
	cerr << "Num NonZero Particles:  " << nonZeroParticle << endl;

	fprintf( stderr, "Complete\n" );
*/
}
