#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>

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
using std::stringstream;

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

#ifndef DO_MPI
#define DO_MPI 0
#else
#include "mpi.h"
#endif


#ifndef DO_VMAG
#define DO_VMAG 0
#endif

void serializeParticles(int start_id, int end_id, int size, double *message, ParticleContainer *advectList);
void serializeParticles(int size, double *mailbox, Particle *particle);
//Arguments:
//  1  2  3  4       5         6       7                8        9
//	Nx Ny Nz ENDTIME PRINT_VTK DO_FTLE ADVECT_PARTICLES STEPSIZE Samples_cubic
//
//
void usage(char* argv[]){
	cout << endl << "Usage for " << argv[0] << "." << endl;
	cout << argv[0] << " accepts the following command line flags along with their default values." << endl << endl;
	cout << "-nx(3)\t\tThe number of points on the domain mesh in the x direction." << endl;	
	cout << "-ny(3)\t\tThe number of points on the domain mesh in the y direction." << endl;
	cout << "-nz(2)\t\tThe number of points on the domain mesh in the z direction." << endl;
	cout << "-t(10)\t\tThe end time of the particle simulation, start time is always 0." << endl;
	cout << "-p(1)\t\tPrint to VTK files or not. Should be 0 (off) or 1 (on)." << endl;
	cout << "-ftle(1)\tPerform an FTLE calculation or not. Should be 0 (off) or 1 (on)." << endl;
	cout << "-a(7)\t\tControls how to seed the mesh with particles and the advection process. Should be one of 0 through 7." << endl;
	cout << "-step(0.001)\tStepsize for the Eulerian algorithm." << endl;
	cout << "-c(10)\t\tNumber of samples in each of x, y and z directions. Requires -a flag to have a value of 7." << endl;
	cout << "-d(0)\t\tWhich data set to run, 0 = small, 1 = Regular, 2 = large." << endl;
	cout << "-u, -usage\tPrints this usage message." << endl;
	cout << endl;
}


int main( int argc, char** argv )
{
	int numProcs, rank;

#if DO_MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#else
	rank = 0;
#endif

	double EndTime = 10;

	int ADVECT_PARTICLES = 7;

	int DO_FTLE = 1;

	int PRINT_VTK = 1;

	double STEPSIZE = 0.001;

	double Samples_cubic = 10;

	int dataSet = 1;

	long int fnx = 3;	
	long int fny = 3;	
	long int fnz = 2;

    double xmin = -118;
    double ymin = -50;
    double zmin = -81.6;
    double xmax = xmin + 168;
    double ymax = ymin + 100;
    double zmax = zmin + 203.6;

	double SampleXmax = xmax;
	double SampleXmin = xmin;
	double SampleYmax = ymax;
	double SampleYmin = ymin;
	double SampleZmax = zmax;
	double SampleZmin = zmin;

	for(int flag = argc; flag > 0; flag--){
		if(!strcmp(argv[flag-1],"-nx")) fnx = atol( argv[flag] );
		if(!strcmp(argv[flag-1],"-ny")) fny = atol( argv[flag] );
		if(!strcmp(argv[flag-1],"-nz")) fnz = atol( argv[flag] );
		if(!strcmp(argv[flag-1],"-t")) EndTime = strtod( argv[flag], NULL );
		if(!strcmp(argv[flag-1],"-p")) PRINT_VTK = atoi( argv[flag] );
		if(!strcmp(argv[flag-1],"-ftle")) DO_FTLE = atoi( argv[flag] );
		if(!strcmp(argv[flag-1],"-a")) ADVECT_PARTICLES = atoi( argv[flag] );
		if(!strcmp(argv[flag-1],"-step")) STEPSIZE = strtod( argv[flag], NULL );
		if(!strcmp(argv[flag-1],"-c")) Samples_cubic = strtod (argv[flag], NULL );
		if(!strcmp(argv[flag-1],"-d")) dataSet = atoi( argv[flag] );
		if(!strcmp(argv[flag-1],"-u") || !strcmp(argv[flag-1],"-usage")){
			usage(argv);
			return 0;
		}
	}

	if( fnx < 2 )
	{
		cerr << "fnx is less than 2, please use a value greater than 2." << endl;
		usage(argv);
		return 1;
	}
	if( fny < 2 )
	{
		cerr << "fny is less than 2, please use a value greater than 2." << endl;
		usage(argv);
		return 1;
	}

	if( fnz < 2 )
	{
		cerr << "fnz is less than 2, please use a value greater than 2." << endl;
		usage(argv);
		return 1;
	}

	if( ADVECT_PARTICLES < 0 || ADVECT_PARTICLES > 7 )
	{
		cerr << "ADVECT_PARTICLES is not valid, should be between 0 and 7 (inclusive)." << endl;
		usage(argv);
		return 1;
	}
	if(PRINT_VTK != 0 && PRINT_VTK != 1){
		cerr << "PRINT_VTK is not valid, should be 0 or 1." << endl;
		usage(argv);
		return 1;
	}

	if( DO_FTLE != 0 && DO_FTLE != 1 )
	{
		cerr << "DO_FTLE is not valid, should be 0 or 1." << endl;
		usage(argv);
		return 1;
	}

	if( STEPSIZE <= 0 )
	{
		cerr << "STEPSIZE should be greater than 0." << endl; 
		usage(argv);
		return 1;
	}
	if(EndTime <= 0 ){
		cerr << "EndTime should be greater than 0." << endl;
		usage(argv);
		return 1;
	}
	
	if( Samples_cubic < 1 )
	{
		cerr << "Samples_cubic should be at least 1." << endl;
		usage(argv);
		return 1;
	}


	cout << "[---------- START ----------]" << endl;

/*
 *Building the problem mesh
 */

	cout << "\tCreating the Problem Mesh" << endl;

// Number of Points in the XYZ dimensions
// 10x10x10 cells
	long int nx, ny, nz;

	if( dataSet == 1 )
	{
		nx = 100;
		ny = 100;
		nz = 100;
	}
	else if( dataSet == 2 )
	{
		nx = 200;
		ny = 200;
		nz = 200;
	}
	else if( dataSet == 3 )
	{
		nx = 300;
		ny = 300;
		nz = 300;
	}
	else if( dataSet == 4 )
	{
		nx = 400;
		ny = 400;
		nz = 400;
	}
	else if( dataSet == 5 )
	{
		nx = 500;
		ny = 500;
		nz = 500;
	}
	else
	{
		nx = 500;
		ny = 500;
		nz = 500;
	}

	long int size = nx*ny*nz*3;
    double* v_field = new double [ size ];
	float *buff = new float [ size ];

// Step Size on Mesh
    double dx = 168 / (nx-1);
    double dy = 100 / (ny-1);
    double dz = 203.6 / (nz-1);

	double BoundBox[6] = { xmin, xmax, ymin, ymax, zmin, zmax };

#if DO_MPI
	if( rank == 0 )
	{
#endif
		cout << "\t\tReading Data Set" << endl;

		long int byteSize = size*sizeof(float);
		
		ifstream is;

		if( dataSet == 1 )
		{
			is.open( "/home/user/Research/DATA/1to5hund/nimrod_100", ios::in | ios::binary );
		}
		else if( dataSet == 2 )
		{
			is.open( "/home/user/Research/DATA/1to5hund/nimrod_200", ios::in | ios::binary );
		}
		else if( dataSet == 3 )
		{
			is.open( "/home/user/Research/DATA/1to5hund/nimrod_300", ios::in | ios::binary );
		}
		else if( dataSet == 4 )
		{
			is.open( "/home/user/Research/DATA/1to5hund/nimrod_400", ios::in | ios::binary );
		}
		else if( dataSet == 5 )
		{
			is.open( "/home/user/Research/DATA/1to5hund/nimrod_500", ios::in | ios::binary );
		}
		else
		{
			is.open( "/home/user/Research/DATA/Fishtank/fishtank", ios::in | ios::binary );
		}

		is.read( reinterpret_cast<char*>(buff), byteSize );

		is.close();
#if DO_MPI
	}

	MPI_Bcast( buff, size, MPI_FLOAT, 0, MPI_COMM_WORLD );

	if( rank == 0 )
	{
#endif
		cout << "\t\tBuilding Velocity Field" << endl;
#if DO_MPI
	}
#endif

#if DO_VMAG
	double magMin = 999999999;
	double magMax = 0.0;
	double magSum = 0.0;
	int numberZeroMag = 0;
#endif

    for( long int z = 0; z < nz; z++ ){
        for( long int y = 0; y < ny; y++ ){
            for( long int x = 0; x < nx; x++ )
            {

                long int id = 3*(x + y*nx + z*nx*ny);

                v_field[ id ]   = buff[id];
                v_field[ id+1 ] = buff[id+1];
                v_field[ id+2 ] = buff[id+2];

#if DO_VMAG
				double mag = buff[id]*buff[id] + buff[id+1]*buff[id+1] + buff[id+2]*buff[id+2];


				if( mag == 0.0 )
				{
					numberZeroMag++;
				}
				else
				{
					magSum += sqrt(mag);
					if( magMin > mag )
					{
						magMin = mag;
					}

					if( magMax < mag )
					{
						magMax = mag;
					}
				}
#endif

            }
        }
    }

#if DO_VMAG
	int numberVs = nx*ny*nz - numberZeroMag;
	double avg = magSum / ( (numberVs==0) ? 1 : numberVs );
	if( avg == magSum ) cerr << "Error" << endl;

	cerr << "Min Mag: " << magMin << endl;
	cerr << "Avg Mag: " << avg << endl;
	cerr << "Max Mag: " << magMax << endl;

	cerr << "AVG: StepSize with Z: " << ((0.8)*((double)(xmax-xmin)/nx) + (0.2)*((double)(zmax-zmin)/nz))/avg << endl;
	cerr << "AVGH: StepSize with Z: " << ((0.8)*((double)(xmax-xmin)/(2*nx)) + (0.2)*((double)(zmax-zmin)/(2*nz)))/avg << endl;
	cerr << "AVG: StepSize no Z  : " << (((double)(xmax-xmin)/nx))/avg << endl;
	cerr << "MIN: StepSize with Z: " << ((0.8)*((double)(xmax-xmin)/nx) + (0.2)*((double)(zmax-zmin)/nz))/magMin << endl;
	cerr << "MIN: StepSize no Z  : " << (((double)(xmax-xmin)/nx))/magMin << endl;
	cerr << "MAX: StepSize with Z: " << ((0.8)*((double)(xmax-xmin)/nx) + (0.2)*((double)(zmax-zmin)/nz))/magMax << endl;
	cerr << "MAX: StepSize no Z  : " << (((double)(xmax-xmin)/nx))/magMax << endl;
 
	exit(0);

#endif

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

cerr << "Init Flows" << endl;
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
			#if DO_MPI
			BuildParticleContainerFull( &FineMesh, advectionList, numParticles, STEPSIZE, rank, numProcs );
			#else
			BuildParticleContainerFull( &FineMesh, advectionList, numParticles, STEPSIZE );
			#endif
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
			double uD[3] = { UberMesh.dx, UberMesh.dy, UberMesh.dz };


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

			double samDx = (SampleXmax-SampleXmin)/(double)numSamples;
			double samDy = (SampleYmax-SampleYmin)/(double)numSamples;
			double samDz = (SampleZmax-SampleZmin)/(double)numSamples;
			double deltaX = 0.01*samDx;
			double deltaY = 0.01*samDy;
			double deltaZ = 0.01*samDz;

			double totalLagrange = 0.0;
			double totalEuler = 0.0;

			double lag_start, lag_end;
			double eul_start, eul_end;

			int totalL = 0;
			int totalE = 0;

cerr << "looping over samples" << endl;
			for( int samplesZ = 0; samplesZ < numSamples; samplesZ++ ){
				for( int samplesY = 0; samplesY < numSamples; samplesY++ ){
					for( int samplesX = 0; samplesX < numSamples; samplesX++ )
					{

						int sampleID = samplesX + samplesY*numSamples + samplesZ*numSamples*numSamples;

						double xp = SampleXmin + samplesX*samDx;
						double yp = SampleYmin + samplesY*samDy;
						double zp = SampleZmin + samplesZ*samDz;
					
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
			#if DO_MPI
			AdvectParticleList( &FineMesh, &UberMesh, &advectList, EndTime, BoundBox, totalL, totalE, rank, numProcs );
			#else	
			AdvectParticleList( &FineMesh, &UberMesh, &advectList, EndTime, BoundBox, totalL, totalE );
			#endif
			GET_TIME (Lagrange_Full_End );
			
			stringstream mpi_cerr;		
	
			mpi_cerr << endl;
			mpi_cerr << rank << ": Total Lagrangian Steps: " << totalL << endl;
			mpi_cerr << rank << ": Total Eulerian Steps:   " << totalE << endl;
			mpi_cerr << rank << ": Average Lagrange PP:    " << (double)totalL / (double)(numParticles) << endl;
			mpi_cerr << rank << ": Average Eulerian PP:    " << (double)totalE / (double)(numParticles) << endl;
			mpi_cerr << rank << ": Percent Lagrange Steps: " << ((double)totalL / (double)( totalL + totalE ) )* 100.0 << endl;

			mpi_cerr << endl;

			cerr << mpi_cerr.rdbuf();

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

					/* We need to pull all particles into global list for everyone to have. ALL Reduce?  */	

					cerr << "Building Advection Lists for Euler only run" << endl;

					Particle* eulerOnly;
					#if DO_MPI
					BuildParticleContainerFull( &FineMesh, eulerOnly, numParticles, STEPSIZE, rank, numProcs );
					#else
					BuildParticleContainerFull( &FineMesh, eulerOnly, numParticles, STEPSIZE );
					#endif
					
					ParticleContainer advectList2( eulerOnly, numParticles );

					cerr << "Advecting euler list" << endl;		

					double Euler_Full_Start;
					double Euler_Full_End;

					GET_TIME( Euler_Full_Start ); 
					AdvectParticleList( &FineMesh, &advectList2, EndTime, BoundBox );
					GET_TIME( Euler_Full_End ); 

					cerr << "building input list for comparison " << endl;

					#if DO_MPI
					Particle* inputPointsL;
					BuildParticleContainerFull( &FineMesh, inputPointsL, numParticles, STEPSIZE, rank, numProcs );
					ParticleContainer advectList3( inputPointsL, numParticles );
					#else
					Particle* inputPoints;
					BuildParticleContainerFull( &FineMesh, inputPoints, numParticles, STEPSIZE );
					#endif

					int totalNumP = nx*ny*nz;

					#if DO_MPI
	
					int start_id = rank*(totalNumP/numProcs);
					int end_id = start_id + (totalNumP/numProcs);
					if( rank == numProcs - 1 )
					{
						end_id += totalNumP - end_id;
					}
	
					MPI_Barrier(MPI_COMM_WORLD);


					Particle* inputPoints = new Particle [totalNumP];
					Particle* outputPoints = new Particle [totalNumP];
					Particle* outputPoints2 = new Particle [totalNumP];

					double *message = new double [4*totalNumP];
					double *mailbox = new double [4*totalNumP];


					memset(mailbox, 0.0, 4*totalNumP*sizeof(double));


					serializeParticles(start_id, end_id, totalNumP, message, &advectList);

					MPI_Allreduce(message, mailbox, totalNumP*4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);						

					serializeParticles(totalNumP, mailbox, outputPoints);


					memset(mailbox, 0.0, 4*totalNumP*sizeof(double));
					serializeParticles(start_id, end_id, totalNumP, message, &advectList2);
					MPI_Allreduce(message, mailbox, totalNumP*4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
					serializeParticles(totalNumP, mailbox, outputPoints2);


					memset(mailbox, 0.0, 4*totalNumP*sizeof(double));
					serializeParticles(start_id, end_id, totalNumP, message, &advectList3);
					MPI_Allreduce(message, mailbox, totalNumP*4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
					serializeParticles(totalNumP, mailbox, inputPoints);

					if( message != NULL)
					{	
						delete [] message;
					}
					if( mailbox != NULL )
					{
						delete [] mailbox;
					}	
					if( inputPointsL != NULL )
					{
						delete [] inputPointsL;
					}

					#else
					
					Particle* outputPoints = advectList.particle;
					Particle* outputPoints2 = advectList2.particle;

					#endif


					
					FTLE1	 = new double [ totalNumP ];
					FTLE2	 = new double [ totalNumP ];
					FTLEDiff = new double [ totalNumP ];

					memset( FTLE1, 0.0, totalNumP*sizeof(double));
					memset( FTLE2, 0.0, totalNumP*sizeof(double));
					memset( FTLEDiff, 0.0, totalNumP*sizeof(double));

					cerr << "Doing FTLE calculation" << endl;

					#if DO_MPI
				
					#pragma omp parallel for schedule(static,1)
					for( int id = start_id; id < end_id; id++ )
					{
						long int ids[3];
						FineMesh.D1to3P( id, ids );
						int x_id = ids[0];
						int y_id = ids[1];
						int z_id = ids[2];

					#else

					#pragma omp parallel for schedule(static,1)
					for( int z_id = 0; z_id < nz; z_id++ ){
						for( int y_id = 0; y_id < ny; y_id++ ){
							for( int x_id = 0; x_id < nx; x_id++ )
							{
								unsigned int id = x_id + y_id * nx + z_id * nx * ny;
					#endif	

								FTLE1[ id ] = compute_FTLE( x_id, y_id, z_id, nx, ny, nz, inputPoints, outputPoints, EndTime ); 
								FTLE2[ id ] = compute_FTLE( x_id, y_id, z_id, nx, ny, nz, inputPoints, outputPoints2, EndTime ); 

								FTLEDiff[ id ] = fabs( FTLE1[id] - FTLE2[id] );

					#if DO_MPI
							}
					#else
							}
						}
					} 
					#endif

					#if DO_MPI
					double *FTLE1_F, *FTLE2_F, *FTLEDiff_F;

					if( rank == 0 )
					{
						FTLE1_F = new double [totalNumP];
						FTLE2_F = new double [totalNumP];
						FTLEDiff_F = new double [totalNumP];

						memset( FTLE1_F, 0.0, totalNumP*sizeof(double));
						memset( FTLE2_F, 0.0, totalNumP*sizeof(double));
						memset( FTLEDiff_F, 0.0, totalNumP*sizeof(double));

					}

					MPI_Reduce( FTLE1, FTLE1_F, totalNumP, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
					MPI_Reduce( FTLE2, FTLE2_F, totalNumP, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
					MPI_Reduce( FTLEDiff, FTLEDiff_F, totalNumP, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

					if( FTLE1 != NULL )
					{
						delete [] FTLE1;	
					}
					if( FTLE2 != NULL )
					{
						delete [] FTLE2;	
					}
					if( FTLEDiff != NULL )
					{
						delete [] FTLEDiff;	
					}
					
					#endif

					#if DO_MPI
					if ( rank == 0 )
					{
						printVtkFTLE( nx, ny, nz, xmin, ymin, zmin, dx, dy, dz, FTLE1_F, FTLE2_F, FTLEDiff_F );
					}
					#else
					printVtkFTLE( nx, ny, nz, xmin, ymin, zmin, dx, dy, dz, FTLE1, FTLE2, FTLEDiff );
					#endif
					

					double Euler_Full = Euler_Full_End - Euler_Full_Start;
					double Lagrange_Full = Lagrange_Full_End - Lagrange_Full_Start;
				
					stringstream mpi_cerr2;

					mpi_cerr2 << rank << ": Timing Full: " << endl;
					mpi_cerr2 << rank << ": Lagrangian: " << Lagrange_Full << endl;
					mpi_cerr2 << rank << ": Eulerian:   " << Euler_Full << endl;
					mpi_cerr2 << rank << ": Lagrange to Euler Full Speedup: " <<  Euler_Full / Lagrange_Full << endl;

					cerr << mpi_cerr2.rdbuf();

				}
			}
		}
	}

	if ( PRINT_VTK != 0 && rank == 0 )
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

#if DO_MPI
	MPI_Finalize();
#endif

}
void serializeParticles(int size, double *mailbox, Particle *particle){
	for(int i = 0; i < size; i++){
		int id = i*4;
		particle[i].x = mailbox[id];
		particle[i].y = mailbox[id+1];
		particle[i].z = mailbox[id+2];
		particle[i].t = mailbox[id+3];
	}
}
void serializeParticles(int start_id, int end_id, int size, double *message, ParticleContainer *advectList){

	for(int i = 0; i < size; i++)
	{
		int id = i*4;
		if(i >= start_id && i < end_id)
		{
			int n_i = i-start_id;
			message[id]   = advectList->particle[n_i].x;
			message[id+1] = advectList->particle[n_i].y;
			message[id+2] = advectList->particle[n_i].z;
			message[id+3] = advectList->particle[n_i].t;
		}
		else
		{
			message[id]   = 0.0;
			message[id+1] = 0.0;
			message[id+2] = 0.0;
			message[id+3] = 0.0;
		}
	}
}

