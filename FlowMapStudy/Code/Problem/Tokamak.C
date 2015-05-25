#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstring>

#include <DomainMesh.h>
#include <Mesh.h>
#include <PrintVTK.h>
#include <AdvectParticles.h>
#include <ParticleContainer.h>
#include <FTLE.h>

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

#define GET_TIME(now){ \
    struct timeval t; \
    gettimeofday(&t, NULL); \
    now = t.tv_sec + t.tv_usec/1000000.0; \
}


//Arguments:
//  1  2  3  4       5         6       7                8        9
//	Nx Ny Nz ENDTIME PRINT_VTK DO_FTLE ADVECT_PARTICLES STEPSIZE Samples_cubic
//


int main( int argc, char** argv )
{

	double EndTime = 10;

	int ADVECT_PARTICLES = 7;

	int DO_FTLE = 1;

	int PRINT_VTK = 1;

	double STEPSIZE = 0.001;

	double Samples_cubic = 10;

	long int fnx = 3;	
	long int fny = 3;	
	long int fnz = 2;

	if( argc >= 4 )
	{
		fnx = atol( argv[1] );	
		fny = atol( argv[2] );	
		fnz = atol( argv[3] );	
	}
	
	if( argc >= 5 )
	{
		EndTime = strtod( argv[4], NULL );
	}

	if( argc >= 6 )
	{
		PRINT_VTK = atoi( argv[5] );	
	}

	if( argc >= 7 )
	{
		DO_FTLE = atoi( argv[6] );	
	}

	if( argc >= 8 )
	{
		ADVECT_PARTICLES = atoi( argv[7] );
	}

	if( argc >= 9 )
	{
		STEPSIZE = strtod( argv[8], NULL );
	}
	
	if( argc >= 10 )
	{
		Samples_cubic = strtod( argv[9], NULL );
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

// Allocate the flows for the flow map
	UberMesh.flowField[0] = new Flow [ num_X_planar_flows ];
	UberMesh.flowField[1] = new Flow [ num_Y_planar_flows ];
	UberMesh.flowField[2] = new Flow [ num_Z_planar_flows ];

// Allocate the accepted flow bool for the flow map
	UberMesh.AcceptableFlow[0] = new bool [ num_X_planar_cells ];
	UberMesh.AcceptableFlow[1] = new bool [ num_Y_planar_cells ];
	UberMesh.AcceptableFlow[2] = new bool [ num_Z_planar_cells ];
// Allocate the max_diff for the flow map
	UberMesh.max_diff[0] = new double [ num_X_planar_cells ]();
	UberMesh.max_diff[1] = new double [ num_Y_planar_cells ]();
	UberMesh.max_diff[2] = new double [ num_Z_planar_cells ]();

//Initisalize flow data
	UberMesh.fillInFlows(); 
	UberMesh.setFlowsCellIDs();

    cout << "\tBuilding Particle List" << endl;

	ParticleContainer flowXParticles( UberMesh.flowField[0], num_X_planar_flows, STEPSIZE );
	ParticleContainer flowYParticles( UberMesh.flowField[1], num_Y_planar_flows, STEPSIZE );
	ParticleContainer flowZParticles( UberMesh.flowField[2], num_Z_planar_flows, STEPSIZE );

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

	cout << "[----------- END -----------]" << endl;


	if( ADVECT_PARTICLES != 0 )
	{
		cout << "Building Advection List" << endl;

		int cbe = 100;

		long int numParticles = 1;
		Particle* advectionList;

		if (ADVECT_PARTICLES == 1 )
			BuildParticleContainerOne( advectionList, STEPSIZE );
		else if (ADVECT_PARTICLES == 2 )
		{
			numParticles = 2;
			BuildParticleContainerTwo( advectionList, STEPSIZE );
		}
		else if (ADVECT_PARTICLES == 3 )
		{
			BuildParticleContainerFull( &FineMesh, advectionList, numParticles, STEPSIZE );
		}
		else if (ADVECT_PARTICLES == 4 )
		{
			BuildParticleContainerFullX2( &FineMesh, advectionList, numParticles, STEPSIZE );
		}
		else if (ADVECT_PARTICLES == 5 )
		{
			numParticles = cbe*cbe*cbe;
			double cube[6] = { -1, 1, -1, 1, -1, 1 };
			BuildParticleContainerCube( advectionList, cube, cbe, STEPSIZE );
		}
		else if( ADVECT_PARTICLES == 6 )
		{

			long int uN[3] = { UberMesh.nx, UberMesh.ny, UberMesh.nz };
			long int uD[3] = { UberMesh.dx, UberMesh.dy, UberMesh.dz };


			cerr << "Building Particles" << endl;

			numParticles = num_X_planar_flows + num_Y_planar_flows + num_Z_planar_flows;
			BuildParticleContainerPlanar( &FineMesh, uN, uD, advectionList, numParticles, STEPSIZE );

			cerr << "Built" << endl;

		}
		else if( ADVECT_PARTICLES == 7 )
		{
			Particle** SampleParticlesI;
			Particle** SampleParticlesE;
			Particle** SampleParticlesL;
			int numSamples = Samples_cubic;

			int SampleSize = numSamples*numSamples*numSamples;

			SampleParticlesI = new Particle* [ SampleSize ];
			SampleParticlesE = new Particle* [ SampleSize ];
			SampleParticlesL = new Particle* [ SampleSize ];

			double FTLE_PD_Avg = 0.0;
			double FTLE_DF_Avg = 0.0;

			int cbe = 3;
			int numParticlesPerCube = 27;

			double samDx = (xmax-xmin)/(double)numSamples;
			double samDy = (ymax-ymin)/(double)numSamples;
			double samDz = (zmax-zmin)/(double)numSamples;
			double deltaX = 0.01*samDx;
			double deltaY = 0.01*samDy;
			double deltaZ = 0.01*samDz;

			double totalLagrange = 0.0;
			double totalEuler = 0.0;

			double lag_start, lag_end;
			double eul_start, eul_end;

			int totalL = 0;
			int totalE = 0;

			for( int samplesZ = 0; samplesZ < numSamples; samplesZ++ ){
				for( int samplesY = 0; samplesY < numSamples; samplesY++ ){
					for( int samplesX = 0; samplesX < numSamples; samplesX++ )
					{

						int sampleID = samplesX + samplesY*numSamples + samplesZ*numSamples*numSamples;

						double xp = xmin + samplesX*samDx;
						double yp = ymin + samplesY*samDy;
						double zp = zmin + samplesZ*samDz;
					
						double cube[6] = { xp-deltaX, xp+deltaX, yp-deltaY, yp+deltaY, zp-deltaZ, zp+deltaZ };

						BuildParticleContainerCube( SampleParticlesI[sampleID], cube, cbe, STEPSIZE );
						BuildParticleContainerCube( SampleParticlesE[sampleID], cube, cbe, STEPSIZE );
						BuildParticleContainerCube( SampleParticlesL[sampleID], cube, cbe, STEPSIZE );

						ParticleContainer advectListE( SampleParticlesE[sampleID], numParticlesPerCube );	
						ParticleContainer advectListL( SampleParticlesL[sampleID], numParticlesPerCube );	

						GET_TIME( eul_start );
						AdvectParticleList( &FineMesh, &advectListE, EndTime, BoundBox );
						GET_TIME( eul_end );

						GET_TIME( lag_start );
						AdvectParticleList( &FineMesh, &UberMesh, &advectListL, EndTime, BoundBox, totalL, totalE );
						GET_TIME( lag_end );

						totalLagrange += lag_end-lag_start;
						totalEuler += eul_end-eul_start;

						Particle* Input   = SampleParticlesI[ sampleID ];
						Particle* OutputE = SampleParticlesE[ sampleID ];
						Particle* OutputL = SampleParticlesL[ sampleID ];

						double FTLE_E, FTLE_L;

						FTLE_E = compute_FTLE( 1, 1, 1, cbe, cbe, cbe, Input, OutputE, EndTime ); 
						FTLE_L = compute_FTLE( 1, 1, 1, cbe, cbe, cbe, Input, OutputL, EndTime ); 

						double avg_ftle = ( FTLE_E + FTLE_L ) / 2.0;

						FTLE_PD_Avg += (( fabs( FTLE_E - FTLE_L ) ) / ( (avg_ftle==0) ? 1 : avg_ftle ) ) * 100;
						FTLE_DF_Avg += fabs( FTLE_E - FTLE_L );

/*
						cerr << "FTLE_E: " << FTLE_E << endl;
						cerr << "FTLE_L: " << FTLE_L << endl;
						cerr << "FTLE_PD: " << (( fabs( FTLE_E - FTLE_L ) ) / ( (avg_ftle==0) ? 1 : avg_ftle ) ) * 100 << endl;
						cerr << "FTLE_DF: " << fabs( FTLE_E - FTLE_L ) << endl;
*/
					}
				}
			}

			FTLE_PD_Avg /= SampleSize;
			FTLE_DF_Avg /= SampleSize;

			cerr << "Num Particles Advected in Total: " << SampleSize*27 << endl;
			cerr << "Num Particles Used in FTLE: " << SampleSize << endl;
			cerr << "Total Lagrangian Steps: " << totalL << endl;
			cerr << "Total Eulerian Steps:   " << totalE << endl;
			cerr << "Average Lagrange PP:    " << (double)totalL / (double)(SampleSize*27) << endl;
			cerr << "Average Eulerian PP:    " << (double)totalE / (double)(SampleSize*27) << endl;
			cerr << "Percent Lagrange Steps: " << ((double)totalL / (double)( totalL + totalE ) )* 100.0 << endl;

			cerr << "Sample FTLE Percent Diff: " << FTLE_PD_Avg << endl;
			cerr << "Sample FTLE Difference:   " << FTLE_DF_Avg << endl;

			cerr << "Total Lag_time:    " << totalLagrange << endl;
			cerr << "Total Eul_time:    " << totalEuler << endl;
			cerr << "Speedup Eul / Lag: " << totalEuler / totalLagrange << endl;


		}

		if( ADVECT_PARTICLES != 7 )
		{

			ParticleContainer advectList( advectionList, numParticles );

			cout << "Advecting Particle List" << endl;
			double Lagrange_Full_Start;
			double Lagrange_Full_End;
			
			int totalL = 0;
			int totalE = 0;

			GET_TIME (Lagrange_Full_Start );
			AdvectParticleList( &FineMesh, &UberMesh, &advectList, EndTime, BoundBox, totalL, totalE );
			GET_TIME (Lagrange_Full_End );
			
			cerr << "Total Lagrangian Steps: " << totalL << endl;
			cerr << "Total Eulerian Steps:   " << totalE << endl;
			cerr << "Average Lagrange PP:    " << (double)totalL / (double)(numParticles) << endl;
			cerr << "Average Eulerian PP:    " << (double)totalE / (double)(numParticles) << endl;
			cerr << "Percent Lagrange Steps: " << ((double)totalL / (double)( totalL + totalE ) )* 100.0 << endl;


			if( DO_FTLE != 0 )
			{

				double *FTLE1;
				double *FTLE2;
				double *FTLEDiff;

				if( ADVECT_PARTICLES != 3 )
				{
					cerr << "Can only do FTLE for full point seeding currently" << endl;
				}
				else
				{
					
					cerr << "Building Advection Lists for Euler only run" << endl;

					Particle* eulerOnly;
					BuildParticleContainerFull( &FineMesh, eulerOnly, numParticles, STEPSIZE );		
					ParticleContainer advectList2( eulerOnly, numParticles );

					cerr << "Advecting euler list" << endl;		

					double Euler_Full_Start;
					double Euler_Full_End;

					GET_TIME( Euler_Full_Start ); 
					AdvectParticleList( &FineMesh, &advectList2, EndTime, BoundBox );
					GET_TIME( Euler_Full_End ); 

					cerr << "building input list for comparison" << endl;

					Particle* inputPoints;
					BuildParticleContainerFull( &FineMesh, inputPoints, numParticles, STEPSIZE );		

					Particle* outputPoints = advectList.particle;
					Particle* outputPoints2 = advectList2.particle;

					FTLE1	 = new double [ nx*ny*nz ];
					FTLE2	 = new double [ nx*ny*nz ];
					FTLEDiff = new double [ nx*ny*nz ];

					cerr << "Doing FTLE calculation" << endl;

					#pragma omp parallel for schedule(static,1)
					for( int z_id = 0; z_id < nz; z_id++ ){
						for( int y_id = 0; y_id < ny; y_id++ ){
							for( int x_id = 0; x_id < nx; x_id++ )
							{
								unsigned int id = x_id + y_id * nx + z_id * nx * ny;

								FTLE1[ id ] = compute_FTLE( x_id, y_id, z_id, nx, ny, nz, inputPoints, outputPoints, EndTime ); 
								FTLE2[ id ] = compute_FTLE( x_id, y_id, z_id, nx, ny, nz, inputPoints, outputPoints2, EndTime ); 

								FTLEDiff[ id ] = fabs( FTLE1[id] - FTLE2[id] );

							}
						}
					} 

					printVtkFTLE( nx, ny, nz, xmin, ymin, zmin, dx, dy, dz, FTLE1, FTLE2, FTLEDiff );

					double Euler_Full = Euler_Full_End - Euler_Full_Start;
					double Lagrange_Full = Lagrange_Full_End - Lagrange_Full_Start;
				

					cerr << "Timing Full: " << endl;
					cerr << "Lagrangian: " << Lagrange_Full << endl;
					cerr << "Eulerian:   " << Euler_Full << endl;
					cerr << "Lagrange to Euler Full Speedup: " <<  Euler_Full / Lagrange_Full << endl;

				}
			}
		}
	}

	if ( PRINT_VTK != 0 )
	{
		cout << "Printing VTK Files" << endl;

		/*VTK File Printer*/
		int umN[3]     = {fnx, fny, fnz};
		int fmN[3]     = {nx, ny, nz};
		double umD[3]  = {fdx, fdy, fdz};
		double fmD[3]  = {dx, dy, dz};
		double orig[3] = { xmin, ymin, zmin };

		VTKFilePrinter printer( umN, fmN, umD, fmD, orig, UberMesh.AcceptableFlow, UberMesh.max_diff, "Tok" );

		printer.printVtkDs();

	}


/*Clean Up of allocatd variables*/


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

