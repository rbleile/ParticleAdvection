#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <fstream>

#include "Flow_Mesh.h"
#include "Mesh.h"
#include "PrintVTK.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::ios;

using std::min;
using std::max;

#include <sys/time.h>

#ifndef DEBUG
#define debug false
#else
#define debug true
#endif


#ifndef STEPSIZE
#define STEPSIZE 0.001
#endif

#ifndef ENDTIME
#define ENDTIME 100.0
#endif


#define GET_TIME(now){ \
    struct timeval t; \
    gettimeofday(&t, NULL); \
    now = t.tv_sec + t.tv_usec/1000000.0; \
}


void AdvectParticleList( Mesh* FineMesh, Flow_Mesh* UberMesh, ParticleList* advectList, double endtime, double* MBB );

void BuildParticleListFullX2( Mesh* mesh, Particle *pl, long int &numP )
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


void BuildParticleListFull( Mesh* mesh, Particle *pl, long int &numP )
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

void BuildParticleListOne( Particle* &pl )
{
	pl = new Particle [1];
	pl[0].x = 1.0;
	pl[0].y = 1.0;
	pl[0].z = 1.0;
	pl[0].setStepSize( STEPSIZE );
}

int main()
{

	cout << "[---------- START ----------]" << endl;

/*
 *Building the problem mesh
 */

	cout << "\tCreating the Problem Mesh" << endl;

// Number of Points in the XYZ dimensions
// 10x10x10 cells
    const long int nx = 300;
    const long int ny = 300;
    const long int nz = 300;

    double* v_field = new double [ 3 * nx * ny * nz ];

// The minimum Corner of the Mesh
    double xmin = -2.44577;
    double ymin = -2.44577;
    double zmin = -1.4629;

// Step Size on Mesh
    double dx = 4.89154 / (nx-1);
    double dy = 4.89154 / (ny-1);
    double dz = 2.92453 / (nz-1);

// The Maximum Corner of the Mesh
    double xmax = -2.44577 + 4.89154;
    double ymax = -2.44577 + 4.89154;
    double zmax = -1.4629  + 2.92453;

	double BoundBox[6] = { xmin, xmax, ymin, ymax, zmin, zmax };

	cout << "\t\tReading Data Set" << endl;

	long int size = nx*ny*nz*3;
	long int byteSize = size*sizeof(float);
	
	ifstream is( "nimrod", ios::in | ios::binary );

	float *buff = new float [ size ];

	is.read( reinterpret_cast<char*>(buff), byteSize );

	is.close();

    cout << "\t\tBuilding Velocity Field" << endl;

    for( long int z = 0; z < nz; z++ ){
        for( long int y = 0; y < ny; y++ ){
            for( long int x = 0; x < nx; x++ )
            {

                long int id = 3*(x + y*nx + z*nx*ny);

                v_field[ id ]   = buff[id];
                v_field[ id+1 ] = buff[id+1];
                v_field[ id+2 ] = buff[id+2];

            }
        }
    }

	if( buff != NULL ) delete [] buff;

    cout << "\tBuilding Fine_Mesh" << endl;
    Mesh FineMesh( nx, ny, nz, dx, dy, dz, xmin, ymin, zmin, v_field );

/*
 * Building a flow map mesh
 */

// Number of points for flow map in the X & Y & Z dims
// 2x2x2 cells
//	good 15
    const long int fnx = 3;
    const long int fny = 3;
    const long int fnz = 3;

// Step Size on Mesh
    double fdx = (xmax-xmin)/(fnx-1);
    double fdy = (ymax-ymin)/(fny-1);
    double fdz = (zmax-zmin)/(fnz-1);

    cout << "\tBuilding Uber_Mesh" << endl;
    Flow_Mesh UberMesh( fnx, fny, fnz, fdx, fdy, fdz, xmin, ymin, zmin, &FineMesh );

    cout << "\t\tAllocating Uber_Mesh Variables" << endl;
// Establish the flows in each planar set
	long int num_X_planar_flows = fnx * ny * nz;
	long int num_Y_planar_flows = fny * nz * nx;
	long int num_Z_planar_flows = fnz * nx * ny;
	long int num_X_planar_cells = (fnx) * (ny-1) * (nz-1);
	long int num_Y_planar_cells = (fny) * (nz-1) * (nx-1);
	long int num_Z_planar_cells = (fnz) * (nx-1) * (ny-1);

/*
cerr << num_X_planar_flows << " " << num_Y_planar_flows << " " << num_Z_planar_flows << endl;
cerr << num_X_planar_cells << " " << num_Y_planar_cells << " " << num_Z_planar_cells << endl;
*/

// Allocate the flows for the flow map
	UberMesh.flowField[0] = new FMFlow [ num_X_planar_flows ];
	UberMesh.flowField[1] = new FMFlow [ num_Y_planar_flows ];
	UberMesh.flowField[2] = new FMFlow [ num_Z_planar_flows ];

// Allocate the accepted flow bool for the flow map
	UberMesh.AcceptableFlow[0] = new bool [ num_X_planar_cells ];
	UberMesh.AcceptableFlow[1] = new bool [ num_Y_planar_cells ];
	UberMesh.AcceptableFlow[2] = new bool [ num_Z_planar_cells ];

//Initisalize flow data
	UberMesh.fillInFlows(); 
	UberMesh.setFlowsCellIDs();

    cout << "\tBuilding Particle List" << endl;

	ParticleList flowXParticles( UberMesh.flowField[0], num_X_planar_flows );
	ParticleList flowYParticles( UberMesh.flowField[1], num_Y_planar_flows );
	ParticleList flowZParticles( UberMesh.flowField[2], num_Z_planar_flows );

	ParticleList* flowParticles[3] = { &flowXParticles, &flowYParticles, &flowZParticles };
	long int numFlows[3] = { num_X_planar_flows, num_Y_planar_flows, num_Z_planar_flows };

	cout << "\tGenerating Precomputed Flow Map" << endl;
	for( int d = 0; d < 3; d++ )
	{
		for( long int i = 0; i < numFlows[d]; i++ )
		{
			FMFlow* flow = &UberMesh.flowField[d][i];
			long int cellID = flow->cellID;
			if( cellID != -1 )
			{
				double cbb[6];
				UberMesh.getCellBounds( cellID, cbb );

				Particle part = flowParticles[d]->particle[i];

				//cerr << "Particle Cell Advection: " << d << " " << i << endl;
				UberMesh.EulerCellAdvection( cellID, 10.0, cbb, part );

				flow->out.setPoint( part.x, part.y, part.z, part.t );
				flow->setFID( cbb );
			}
		}
	}

	cout << "\tComputing Acceptable flows" << endl;
	long int num_accepted = UberMesh.computeAllAcceptableFlows();

	cout << num_X_planar_cells << " " <<  num_Y_planar_cells << " " <<  num_Z_planar_cells << endl;

	long int total = num_X_planar_cells + num_Y_planar_cells + num_Z_planar_cells; 

	cout << "Statistics" << endl;
	cout << "Number of Computable Faces: " << num_accepted << endl;
	cout << "Number of total  Faces:     " << total << endl;

	cout << "Percentage of computable faces: " << ((double)num_accepted / (double)total ) * 100  << "%" << endl;

	cout << "[----------- END -----------]" << endl;

	cout << "Printing VTK Files" << endl;

/*VTK File Printer*/
int umN[3]     = {fnx, fny, fnz};
int fmN[3]     = {nx, ny, nz};
double umD[3]  = {fdx, fdy, fdz};
double fmD[3]  = {dx, dy, dz};
double orig[3] = { xmin, ymin, zmin };

VTKFilePrinter printer( umN, fmN, umD, fmD, orig, UberMesh.AcceptableFlow, "Tok" );

printer.printVtkDs();

	cout << "Building Advection List" << endl;

	long int numParticles = 1;
	Particle* advectionList;

	BuildParticleListOne( advectionList );

//	BuildParticleListFull( &FineMesh, advectionList, numParticles );
//	BuildParticleListFullX2( &FineMesh, advectionList, numParticles );
//

	ParticleList advectList( advectionList, numParticles );

	AdvectParticleList( &FineMesh, &UberMesh, &advectList, ENDTIME, BoundBox );

	if( UberMesh.flowField[0] != NULL )
	{
		delete [] UberMesh.flowField[0];
	}
	if( UberMesh.flowField[1] != NULL )
	{
		delete [] UberMesh.flowField[1];
	}
	if( UberMesh.flowField[2] != NULL )
	{
		delete [] UberMesh.flowField[2];
	}

	if( UberMesh.AcceptableFlow[0] != NULL )
	{
		delete [] UberMesh.AcceptableFlow[0];
	}
	if( UberMesh.AcceptableFlow[1] != NULL )
	{
		delete [] UberMesh.AcceptableFlow[1];
	}
	if( UberMesh.AcceptableFlow[2] != NULL )
	{
		delete [] UberMesh.AcceptableFlow[2];
	}

}

long int GetParticleCellID( Mesh* FineMesh, Flow_Mesh *UberMesh, Particle* part )
{
	double vel[3];
	FineMesh->getVelocity( *part, vel );

	double xmax = FineMesh->x0 + FineMesh->dx*(FineMesh->nx-1);
	double ymax = FineMesh->y0 + FineMesh->dy*(FineMesh->ny-1);
	double zmax = FineMesh->z0 + FineMesh->dz*(FineMesh->nz-1);

	FMPoint np( part->x + vel[0]*STEPSIZE, part->y + vel[1]*STEPSIZE, part->z + vel[2]*STEPSIZE );

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

void ComputeFaceID( long int cellID, Flow_Mesh *UberMesh, Particle* part )
{
	double bb[6];
	UberMesh->getCellBounds( cellID, bb );
	part->setFID( bb );
}

bool AcceptableFlow( long int cellID, Flow_Mesh *UberMesh, Mesh* FineMesh, Particle* part )
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

int AdvectParticleOnFlow( int cellID, Mesh *FineMesh, Flow_Mesh *UberMesh, Particle* part )
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

		FMFlow* flows = UberMesh->flowField[d];

		FMPoint p0 = flows[ flowIDs[0] ].in;
		FMPoint p1 = flows[ flowIDs[1] ].in;
		FMPoint p2 = flows[ flowIDs[2] ].in;
		FMPoint p3 = flows[ flowIDs[3] ].in;

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

void AdvectParticleList( Mesh *FineMesh, Flow_Mesh *UberMesh, ParticleList* advectList, double endtime, double* MBB )
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
