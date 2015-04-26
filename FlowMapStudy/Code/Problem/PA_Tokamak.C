#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstring>

#include "DomainMesh.h"
#include "Mesh.h"
#include "PrintVTK.h"
#include "AdvectParticles.h"
#include "ParticleContainer.h"

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


int main( int argc, char** argv )
{

	long int fnx, fny, fnz;	

	if( argc < 4 )
	{
		fnx = 3;	
		fny = 3;	
		fnz = 3;	
	}
	else
	{
		fnx = atol( argv[1] );	
		fny = atol( argv[2] );	
		fnz = atol( argv[3] );	
	}

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
	
	ifstream is( "/home/user/Research/DATA/nimrod", ios::in | ios::binary );

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
//    const long int fnx = 15;
//    const long int fny = 15;
//    const long int fnz = 15;

// Step Size on Mesh
    double fdx = (xmax-xmin)/(fnx-1);
    double fdy = (ymax-ymin)/(fny-1);
    double fdz = (zmax-zmin)/(fnz-1);

    cout << "\tBuilding Uber_Mesh" << endl;
    DomainMesh UberMesh( fnx, fny, fnz, fdx, fdy, fdz, xmin, ymin, zmin, &FineMesh );

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
	UberMesh.flowField[0] = new Flow [ num_X_planar_flows ];
	UberMesh.flowField[1] = new Flow [ num_Y_planar_flows ];
	UberMesh.flowField[2] = new Flow [ num_Z_planar_flows ];

// Allocate the accepted flow bool for the flow map
	UberMesh.AcceptableFlow[0] = new bool [ num_X_planar_cells ];
	UberMesh.AcceptableFlow[1] = new bool [ num_Y_planar_cells ];
	UberMesh.AcceptableFlow[2] = new bool [ num_Z_planar_cells ];

//Initisalize flow data
	UberMesh.fillInFlows(); 
	UberMesh.setFlowsCellIDs();

    cout << "\tBuilding Particle List" << endl;

	ParticleContainer flowXParticles( UberMesh.flowField[0], num_X_planar_flows );
	ParticleContainer flowYParticles( UberMesh.flowField[1], num_Y_planar_flows );
	ParticleContainer flowZParticles( UberMesh.flowField[2], num_Z_planar_flows );

	ParticleContainer* flowParticles[3] = { &flowXParticles, &flowYParticles, &flowZParticles };
	long int numFlows[3] = { num_X_planar_flows, num_Y_planar_flows, num_Z_planar_flows };

	cout << "\tGenerating Precomputed Flow Map" << endl;
	double start,end;

	GET_TIME( start );

	for( int d = 0; d < 3; d++ )
	{
		#pragma omp parallel for
		for( long int i = 0; i < numFlows[d]; i++ )
		{
			Flow* flow = &UberMesh.flowField[d][i];
			long int cellID = flow->cellID;
			if( cellID != -1 )
			{
				double cbb[6];
				UberMesh.getCellBounds( cellID, cbb );

				Particle part = flowParticles[d]->particle[i];

				//cerr << "Particle Cell Advection: " << d << " " << i << endl;
				UberMesh.EulerCellAdvection( cellID, 100.0, cbb, part );

				flow->out.setPoint( part.x, part.y, part.z, part.t );
				flow->setFID( cbb );
			}
		}
	}

	GET_TIME( end );

	cout << "\t\tTime To Compute FlowMap: " << end-start << " (seconds)" << endl;

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

#if 0
	cout << "Building Advection List" << endl;

	long int numParticles = 1;
	Particle* advectionList;

	BuildParticleContainerOne( advectionList );
//	BuildParticleContainerFull( &FineMesh, advectionList, numParticles );
//	BuildParticleContainerFullX2( &FineMesh, advectionList, numParticles );

	ParticleContainer advectList( advectionList, numParticles );

	AdvectParticleList( &FineMesh, &UberMesh, &advectList, ENDTIME, BoundBox );

#endif

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

