#include <FTLE.h>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cmath>

using namespace Eigen;
using namespace std;

// Distance Between two points
double dist( Particle p1, Particle p2 )
{
	return sqrt( (p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y) + (p1.z-p2.z)*(p1.z-p2.z) );
}

// returns the ftle for a given 3d index
double compute_FTLE( int x_id, int y_id, int z_id, const long int nx, const long int ny, const long int nz, Particle* inputs, Particle* outputs, double endTime )
{

	// If we are not on the boundary of the problem
	if( (x_id)*(x_id+1-nx) < 0 && (y_id)*(y_id+1-ny) < 0 && (z_id)*(z_id+1-nz) < 0 )
	{

		// Id's of the 6 nearest neighbors points
		long int id10 = (x_id-1) + y_id * nx + z_id * nx * ny;
		long int id11 = (x_id+1) + y_id * nx + z_id * nx * ny;
		long int id20 = x_id + (y_id-1) * nx + z_id * nx * ny;
		long int id21 = x_id + (y_id+1) * nx + z_id * nx * ny;
		long int id30 = x_id + y_id * nx + (z_id-1) * nx * ny;
		long int id31 = x_id + y_id * nx + (z_id+1) * nx * ny;

		// Matrix Values for eigen value check
		double A1 = ( dist( outputs[id11], outputs[id10] )  ) / ( dist ( inputs[id11], inputs[id10] ) );
		double A2 = ( dist( outputs[id21], outputs[id20] )  ) / ( dist ( inputs[id21], inputs[id20] ) );
		double A3 = ( dist( outputs[id31], outputs[id30] )  ) / ( dist ( inputs[id31], inputs[id30] ) );

		// Building Matrix A'*A
		Matrix3d A;
		A << A1, A2, A3,
			 A1, A2, A3,
			 A1, A2, A3;

		Matrix3d At;
		At << A1, A1, A1,
			  A2, A2, A2,
			  A3, A3, A3;

		Matrix3d B;
		B = At*A;

		EigenSolver<MatrixXd> es(B);

		// Get the eigen values from B
		Vector3d eigenValues;
		eigenValues <<	es.eigenvalues()[0].real(),
						es.eigenvalues()[1].real(),
						es.eigenvalues()[2].real();

		// Compute the maximum eigenValue

		double delta = eigenValues.maxCoeff();

		if( delta == 0.0 ){ return 0.0; }

//Given endTime is in units where start time is 0. else do endTime-startTime		
		// return value for ftle computation
		return ( log( delta ) / (2*endTime) );
	}
	
	return 0.0;

}


