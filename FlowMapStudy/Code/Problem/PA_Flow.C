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

using std::min;
using std::max;

#include <sys/time.h>

#ifndef DEBUG
#define debug false
#else
#define debug true
#endif

#define GET_TIME(now){ \
    struct timeval t; \
    gettimeofday(&t, NULL); \
    now = t.tv_sec + t.tv_usec/1000000.0; \
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
    const int nx = 111;
    const int ny = 111;
    const int nz = 111;

    double* v_field = new double [ 3 * nx * ny * nz ];

// The minimum Corner of the Mesh
    double xmin = -10.0;
    double ymin = -10.0;
    double zmin = -10.0;

// The Maximum Corner of the Mesh
    double xmax = 10.0;
    double ymax = 10.0;
    double zmax = 10.0;

// Step Size on Mesh
    double dx = (xmax-xmin)/(nx-1);
    double dy = (ymax-ymin)/(ny-1);
    double dz = (zmax-zmin)/(nz-1);

	double BoundBox[6] = { xmin, xmax, ymin, ymax, zmin, zmax };

    cout << "\t\tBuilding Velocity Field" << endl;

    for( int z = 0; z < nz; z++ ){
        for( int y = 0; y < ny; y++ ){
            for( int x = 0; x < nx; x++ )
            {

				double xp = xmin + x*dx;
				double yp = ymin + y*dy;
				double zp = zmin + z*dz;

                double vx = 1;
                double vy = 0;
                double vz = 0;

				double mag = sqrt( (vx*vx) + (vy*vy) + (vz*vz) );

				vx /= mag;
				vy /= mag;
				vz /= mag;
	
                int id = 3*(x + y*nx + z*nx*ny);
                
                v_field[  id  ] = vx;
                v_field[ id+1 ] = vy;
                v_field[ id+2 ] = vz;

				//cerr << "id: " << id << "\t <" << vx << ", " << vy << ", " << vz << " >" << endl;

            }
        }
    }

    cout << "\tBuilding Fine_Mesh" << endl;
    Mesh FineMesh( nx, ny, nz, dx, dy, dz, xmin, ymin, zmin, v_field );

/*
 * Building a flow map mesh
 */

// Number of points for flow map in the X & Y & Z dims
// 2x2x2 cells
    const int fnx = 3;
    const int fny = 3;
    const int fnz = 3;

// Step Size on Mesh
    double fdx = (xmax-xmin)/(fnx-1);
    double fdy = (ymax-ymin)/(fny-1);
    double fdz = (zmax-zmin)/(fnz-1);

    cout << "\tBuilding Uber_Mesh" << endl;
    Flow_Mesh UberMesh( fnx, fny, fnz, fdx, fdy, fdz, xmin, ymin, zmin, &FineMesh );

    cout << "\t\tAllocating Uber_Mesh Variables" << endl;
// Establish the flows in each planar set
	int num_X_planar_flows = fnx * ny * nz;
	int num_Y_planar_flows = fny * nx * nz;
	int num_Z_planar_flows = fnz * nx * ny;
	int num_X_planar_cells = (fnx) * (ny-1) * (nz-1);
	int num_Y_planar_cells = (fny) * (nx-1) * (nz-1);
	int num_Z_planar_cells = (fnz) * (nx-1) * (ny-1);

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
	int numFlows[3] = { num_X_planar_flows, num_Y_planar_flows, num_Z_planar_flows };

	cout << "\tGenerating Precomputed Flow Map" << endl;
	for( int d = 0; d < 3; d++ )
	{
		for( int i = 0; i < numFlows[d]; i++ )
		{
			FMFlow* flow = &UberMesh.flowField[d][i];
			int cellID = flow->cellID;
			if( cellID != -1 )
			{
				double cbb[6];
				UberMesh.getCellBounds( cellID, cbb );

				Particle part = flowParticles[d]->particle[i];

				//cerr << "Particle Cell Advection: " << d << " " << i << endl;
				UberMesh.EulerCellAdvection( cellID, 1000.0, cbb, part );

				flow->out.setPoint( part.x, part.y, part.z, part.t );
				flow->setFID( cbb );
			}
		}
	}

	cout << "\tComputing Acceptable flows" << endl;
	int num_accepted = UberMesh.computeAllAcceptableFlows();

	//cerr << num_X_planar_cells << " " <<  num_Y_planar_cells << " " <<  num_Z_planar_cells << endl;

	int total = num_X_planar_cells + num_Y_planar_cells + num_Z_planar_cells; 

	cout << "Statistics" << endl;
	cout << "Number of Computable Faces: " << num_accepted << endl;
	cout << "Number of total  Faces:     " << total << endl;

	cout << "Percentage of computable faces: " << ((double)num_accepted / (double)total ) * 100  << "%" << endl;

	cout << "[----------- END -----------]" << endl;


/*VTK File Printer*/
int umN[3]     = {fnx, fny, fnz};
int fmN[3]     = {nx, ny, nz};
double umD[3]  = {fdx, fdy, fdz};
double fmD[3]  = {dx, dy, dz};
double orig[3] = { xmin, ymin, zmin };

VTKFilePrinter printer( umN, fmN, umD, fmD, orig, UberMesh.AcceptableFlow, "Cus" );

printer.printVtkDs();


/*
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
*/
//here
#if 0
	
	/* Start Advection Using Flow Map */
	/* Advection Pattern - Face Centers*/

	double advectionTime = 100.0;

	for( int l = 0; l < 3; l++ )
	{
		int npx = nx;
		int npy = ny;
		int npz = nz;

		double px0 = xmin;
		double py0 = ymin;
		double pz0 = zmin;

		if( l == 0 )
		{
			npy--;
			npz--;

			py0 += dy/2.0;
			pz0 += dz/2.0;

		}

		Particle* particles = new Particle [ npy*npx*npz ];

		for( int k = 0; k < npz; k++ )  // Loop over Z
		{
			for( int j = 0; j < npy; j++ ) // Loop over Y
			{
				for( int i = 0; i < npx; i++ ) // Loop over Z
				{
					int index = i + j*npx + k*npx*npy;

					particles[index].x = px0 + i*dx;
					particles[index].y = py0 + j*dy;
					particles[index].z = pz0 + k*dz;

				}
			}
		}

		ParticleList pl( particles, npx*npy*npz );

		trackParticlesFlow( &mesh, &pl, BoundBox, advectionTime );

	}
#endif
}
