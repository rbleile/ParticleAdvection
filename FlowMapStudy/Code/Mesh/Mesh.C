//Local Includes
#include "Mesh.h"
#include "Flow.h"

//Library Includes
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>

//Name Spaces
using std::cerr;
using std::cout;
using std::endl;
using std::ofstream;

using std::min;
using std::max;

/*System Based Timing*/
#include <sys/time.h>

#define GET_TIME(now){ \
    struct timeval t; \
    gettimeofday(&t, NULL); \
    now = t.tv_sec + t.tv_usec/1000000.0; \
}

/* Define Debug Mode <on/off> ( if debug is in the compile line ) */
#ifndef DEBUG
#define debug false
#else
#define debug true
#endif

//Function for checking floating point equvalence to epsilon error
inline bool almostEqual( float ip1, float ip2 )
{
	float p1 = fabs( ip1 );
	float p2 = fabs( ip2 );
	if( p1 < 1 || p2 < 1 )
	{
		return ( fabs( ip1-ip2 ) <= epsilon );
	}

	return ( ( fabs((ip1-ip2)) <= (( ( p1 < p2) ? p2 : p1 )*epsilon) ) ? true : false);
} 
inline bool almostEqual( double ip1, double ip2 )
{
	double p1 = fabs( ip1 );
	double p2 = fabs( ip2 );
	if( p1 < 1 || p2 < 1 )
	{
		return ( fabs( ip1-ip2 ) <= epsilon );
	}

	return ( ( fabs((ip1-ip2)) <= (( ( p1 < p2) ? p2 : p1 )*epsilon) ) ? true : false);
} 

/***************************************************************************************************
-
-	getVelocity
-		return the velocity at a point in the mesh
-			Trilinear Interpolation of velocity field on the mesh
-
***************************************************************************************************/

void Mesh::getVelocity( Particle point, double *outVelocity ){

	long int pid000[3];
	double px = point.x;
	double py = point.y;
	double pz = point.z;

	getLogicalCellID( px, py, pz, pid000 );

	if( pid000[0] == nx-1 )
	{
	
		if( almostEqual( dx*(nx-1) + x0, px) )
		{
			pid000[0]--;
		}
	}

	if( pid000[1] == ny-1 )
	{
		if( almostEqual( dy*(ny-1) + y0, py) )
		{
			pid000[1]--;
		}
	}

	if( pid000[2] == nz-1 )
	{
		if( almostEqual( dz*(nz-1) + z0, pz) )
		{
			pid000[2]--;
		}
	}

	if( pid000[0] > nx-2 || pid000[1] > ny-2 || pid000[2] > nz-2 )
	{

/*
		cerr << "getVelocity: Error: pid out of bounds" << endl;
		cerr << nx-2 << " " << ny-2 << " " << nz-2 << endl;
		cerr << "Particle: " << px << " " << py << " " << pz << endl;
		cerr << "PID:      " << pid000[0] << " " << pid000[1] << " " << pid000[2] << endl;
		cerr << "Why!!" << endl;
*/

		outVelocity[0] = 0;
		outVelocity[1] = 0;
		outVelocity[2] = 0;
		return;
	}

	double x00 = dx*pid000[0] + x0;
	double y00 = dy*pid000[1] + y0;
	double z00 = dz*pid000[2] + z0;

	double x11 = dx*(pid000[0]+1) + x0;
	double y11 = dy*(pid000[1]+1) + y0;
	double z11 = dz*(pid000[2]+1) + z0;

	double xd = ( ( px - x00 ) / ( x11 - x00 ) );
	double yd = ( ( py - y00 ) / ( y11 - y00 ) );
	double zd = ( ( pz - z00 ) / ( z11 - z00 ) );

	long int id000 = pid000[0] + pid000[1]*( nx  ) + pid000[2]*( nx )*( ny );
	long int id001 = 3*(id000 + 1);
	long int id010 = 3*(id000 + (nx));
	long int id011 = 3*(id000 + (nx) + 1);

	long int id100 = (id000 + (nx)*(ny));
	long int id101 = 3*(id100 + 1);
	long int id110 = 3*(id100 + (nx));
	long int id111 = 3*(id100 + (nx) + 1);
	
	id000 *= 3;
	id100 *= 3;

	double c00x = ( ( V[id000] )   * (1 - xd) ) + ( ( V[id001])   * xd );
	double c00y = ( ( V[1+id000] ) * (1 - xd) ) + ( ( V[1+id001]) * xd );
	double c00z = ( ( V[2+id000] ) * (1 - xd) ) + ( ( V[2+id001]) * xd );

	double c10x = ( ( V[id010] )   * (1 - xd) ) + ( ( V[id011])   * xd );
	double c10y = ( ( V[1+id010] ) * (1 - xd) ) + ( ( V[1+id011]) * xd );
	double c10z = ( ( V[2+id010] ) * (1 - xd) ) + ( ( V[2+id011]) * xd );

	double c01x = ( ( V[id100] )   * (1 - xd) ) + ( ( V[id101])   * xd );
	double c01y = ( ( V[1+id100] ) * (1 - xd) ) + ( ( V[1+id101]) * xd );
	double c01z = ( ( V[2+id100] ) * (1 - xd) ) + ( ( V[2+id101]) * xd );

	double c11x = ( ( V[id110] )   * (1 - xd) ) + ( ( V[id111])   * xd );
	double c11y = ( ( V[1+id110] ) * (1 - xd) ) + ( ( V[1+id111]) * xd );
	double c11z = ( ( V[2+id110] ) * (1 - xd) ) + ( ( V[2+id111]) * xd );

	double c0x = c00x * ( 1 - yd ) + c10x * yd;
	double c0y = c00y * ( 1 - yd ) + c10y * yd;
	double c0z = c00z * ( 1 - yd ) + c10z * yd;

	double c1x = c01x * ( 1 - yd ) + c11x * yd;
	double c1y = c01y * ( 1 - yd ) + c11y * yd;
	double c1z = c01z * ( 1 - yd ) + c11z * yd;

	outVelocity[0] = c0x * ( 1 - zd ) + c1x * zd;    
	outVelocity[1] = c0y * ( 1 - zd ) + c1y * zd;    
	outVelocity[2] = c0z * ( 1 - zd ) + c1z * zd;    

}

void Mesh::getVelocity( Point point, double *outVelocity )
{

	long int pid000[3];
	double px = point.x;
	double py = point.y;
	double pz = point.z;

	getLogicalCellID( px, py, pz, pid000 );

	if( pid000[0] == nx-1 )
	{
	
		if( almostEqual( dx*(nx-1) + x0, px) )
		{
			pid000[0]--;
		}
	}

	if( pid000[1] == ny-1 )
	{
		if( almostEqual( dy*(ny-1) + y0, py) )
		{
			pid000[1]--;
		}
	}

	if( pid000[2] == nz-1 )
	{
		if( almostEqual( dz*(nz-1) + z0, pz) )
		{
			pid000[2]--;
		}
	}

	if( pid000[0] > nx-2 || pid000[1] > ny-2 || pid000[2] > nz-2 )
	{
/*
		cerr << "getVelocity: Error: pid out of bounds" << endl;
		cerr << nx-2 << " " << ny-2 << " " << nz-2 << endl;
		cerr << "Particle: " << px << " " << py << " " << pz << endl;
		cerr << "PID:      " << pid000[0] << " " << pid000[1] << " " << pid000[2] << endl;
		cerr << "Why!!" << endl;
*/
		outVelocity[0] = 0;
		outVelocity[1] = 0;
		outVelocity[2] = 0;
		return;
	}

	double x00 = dx*pid000[0] + x0;
	double y00 = dy*pid000[1] + y0;
	double z00 = dz*pid000[2] + z0;

	double x11 = dx*(pid000[0]+1) + x0;
	double y11 = dy*(pid000[1]+1) + y0;
	double z11 = dz*(pid000[2]+1) + z0;

	double xd = ( ( px - x00 ) / ( x11 - x00 ) );
	double yd = ( ( py - y00 ) / ( y11 - y00 ) );
	double zd = ( ( pz - z00 ) / ( z11 - z00 ) );

	long int id000 = pid000[0] + pid000[1]*( nx  ) + pid000[2]*( nx )*( ny );
	long int id001 = 3*(id000 + 1);
	long int id010 = 3*(id000 + (nx));
	long int id011 = 3*(id000 + (nx) + 1);

	long int id100 = (id000 + (nx)*(ny));
	long int id101 = 3*(id100 + 1);
	long int id110 = 3*(id100 + (nx));
	long int id111 = 3*(id100 + (nx) + 1);
	
	id000 *= 3;
	id100 *= 3;

	double c00x = ( ( V[id000] )   * (1 - xd) ) + ( ( V[id001])   * xd );
	double c00y = ( ( V[1+id000] ) * (1 - xd) ) + ( ( V[1+id001]) * xd );
	double c00z = ( ( V[2+id000] ) * (1 - xd) ) + ( ( V[2+id001]) * xd );

	double c10x = ( ( V[id010] )   * (1 - xd) ) + ( ( V[id011])   * xd );
	double c10y = ( ( V[1+id010] ) * (1 - xd) ) + ( ( V[1+id011]) * xd );
	double c10z = ( ( V[2+id010] ) * (1 - xd) ) + ( ( V[2+id011]) * xd );

	double c01x = ( ( V[id100] )   * (1 - xd) ) + ( ( V[id101])   * xd );
	double c01y = ( ( V[1+id100] ) * (1 - xd) ) + ( ( V[1+id101]) * xd );
	double c01z = ( ( V[2+id100] ) * (1 - xd) ) + ( ( V[2+id101]) * xd );

	double c11x = ( ( V[id110] )   * (1 - xd) ) + ( ( V[id111])   * xd );
	double c11y = ( ( V[1+id110] ) * (1 - xd) ) + ( ( V[1+id111]) * xd );
	double c11z = ( ( V[2+id110] ) * (1 - xd) ) + ( ( V[2+id111]) * xd );

	double c0x = c00x * ( 1 - yd ) + c10x * yd;
	double c0y = c00y * ( 1 - yd ) + c10y * yd;
	double c0z = c00z * ( 1 - yd ) + c10z * yd;

	double c1x = c01x * ( 1 - yd ) + c11x * yd;
	double c1y = c01y * ( 1 - yd ) + c11y * yd;
	double c1z = c01z * ( 1 - yd ) + c11z * yd;

	outVelocity[0] = c0x * ( 1 - zd ) + c1x * zd;    
	outVelocity[1] = c0y * ( 1 - zd ) + c1y * zd;    
	outVelocity[2] = c0z * ( 1 - zd ) + c1z * zd;    

}

/***************************************************************************************************
-
-	stepBackToBoundary
-		Step a particle back to the boundry of a cell to test for flows
-
***************************************************************************************************/
inline void stepBackToBoundry( Particle &particle, double* bbox, double *vel, double dt )
{
	double dtx = 0; double dty = 0; double dtz = 0;

	if( particle.x < bbox[0] )
	{
		double xp = bbox[0] - particle.x;
		double vdt = particle.x - ( particle.x - vel[0]*dt);
		if( vdt != 0.0 ) 
		dtx = (xp / vdt)*dt;
	}
	else if( particle.x > bbox[1] )
	{
		double xp = bbox[1] - particle.x;
		double vdt = particle.x - ( particle.x - vel[0]*dt);
		if( vdt != 0.0 ) 
		dtx = (xp / vdt )*dt;
	}

	if( particle.y < bbox[2] )
	{
		double yp = bbox[2] - particle.y;
		double vdt = particle.y - ( particle.y - vel[1]*dt);
		if( vdt != 0.0 ) 
		dty = (yp / vdt)*dt;
	}
	else if( particle.y > bbox[3] )
	{
		double yp = bbox[3] - particle.y;
		double vdt = particle.y - ( particle.y - vel[1]*dt);
		if( vdt != 0.0 ) 
		dty = (yp / vdt)*dt;
	}

	if( particle.z < bbox[4] )
	{
		double zp = bbox[4] - particle.z;
		double vdt = particle.z - ( particle.z - vel[2]*dt);
		if( vdt != 0.0 ) 
		dtz = ( zp / vdt )*dt;
	}
	else if( particle.z > bbox[5] )
	{
		double zp = bbox[5] - particle.z;
		double vdt = particle.z - ( particle.z - vel[2]*dt);
		if( vdt != 0.0 ) 
		dtz = ( zp / vdt )*dt;
	}
	

	double adtx = fabs( dtx );
	double adty = fabs( dty );
	double adtz = fabs( dtz );
	double ndt = 0;

	if( adtx >= adty && adtx >= adtz ) // x greater than y and z
	{
		ndt = dtx;
	}
	else if( adty >= adtx && adty >= adtz ) // y greater than x and z
	{
		ndt = dty;
	}
	else // z greater than x and y
	{
		ndt = dtz;
	}

	particle.x += vel[0]*ndt;
	particle.y += vel[1]*ndt;
	particle.z += vel[2]*ndt;
	particle.t -= fabs( ndt );


	if( almostEqual( particle.x, bbox[0] ))
	{
		particle.x = bbox[0];
	}
	else if( almostEqual( particle.x, bbox[1] ))
	{
		particle.x = bbox[1];
	}

	if( almostEqual( particle.y, bbox[2] ))
	{
		particle.y = bbox[2];
	}
	else if( almostEqual( particle.y, bbox[3] ))
	{
		particle.y = bbox[3];
	}

	if( almostEqual( particle.z, bbox[4] ))
	{
		particle.z = bbox[4];
	}
	else if( almostEqual( particle.z, bbox[5] ))
	{
		particle.z = bbox[5];
	}
}

/***************************************************************************************************
-
-	Check if the particle is inside a given bounding box
-
***************************************************************************************************/
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

int Mesh::onBoundary( Particle particle, double* bbox )
{

    if( almostEqual( particle.x, bbox[0]) || almostEqual( particle.x, bbox[1] ) ||
        almostEqual( particle.y, bbox[2]) || almostEqual( particle.y, bbox[3] ) ||
        almostEqual( particle.z, bbox[4]) || almostEqual( particle.z, bbox[5]) )
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

inline int onX( Particle particle, double* bbox )
{
	if( (almostEqual( particle.x, bbox[0] ) || almostEqual( particle.x, bbox[1] )) &&
        particle.y >= bbox[2] && particle.y <= bbox[3] &&
        particle.z >= bbox[4] && particle.z <= bbox[5] )
		return 1;
	return 0;
}
inline int onY( Particle particle, double* bbox )
{
	if( (almostEqual( particle.y, bbox[2] ) || almostEqual( particle.y, bbox[3] )) &&
        particle.x >= bbox[0] && particle.x <= bbox[1] &&
        particle.z >= bbox[4] && particle.z <= bbox[5] )
		return 1;
	return 0;
}
inline int onZ( Particle particle, double* bbox )
{
	if( (almostEqual( particle.z, bbox[4] ) || almostEqual( particle.z, bbox[5] )) &&
        particle.y >= bbox[2] && particle.y <= bbox[3] &&
        particle.x >= bbox[0] && particle.x <= bbox[1] )
		return 1;
	return 0;
}

/***************************************************************************************************
-
-	CheckStep
-		Check the boundaries of our step and determine if we need to kill, or step back to bounds
-		return codes ( 0 - kill, 1 - keep going, 2 - reevaluate flows )
-
***************************************************************************************************/
inline int checkStep( Particle &particle, double endTime, double* bbox, double* mBB, double* vel, double dt, int toPrint )
{

	if( toPrint )
	{
		cerr << "Out Pos: " << particle.x << " " << particle.y << " " << particle.z << endl; 
		cerr << "Out Vel: " << vel[0] << " " << vel[1] << " " << vel[2] << endl;
		cerr << "Bbox:    " << bbox[0] << " " << bbox[1] << " \t " << bbox[2] << " " << bbox[3] << " \t " << bbox[4] << " " << bbox[5] << endl;
	}

	if( !checkBoundary(particle, bbox) ){ 

		int inMesh = checkBoundary(particle, mBB); // Check if particle has left mesh entirely

		stepBackToBoundry( particle, bbox, vel, dt );

		if( toPrint )
		{
			cerr << "Out Pos: " << particle.x << " " << particle.y << " " << particle.z << endl; 
			cerr << "Out Vel: " << vel[0] << " " << vel[1] << " " << vel[2] << endl;
			cerr << "Bbox:    " << bbox[0] << " " << bbox[1] << " \t " << bbox[2] << " " << bbox[3] << " \t " << bbox[4] << " " << bbox[5] << endl;
		}

		if( !inMesh ) return 0;
		else          return 2;

	}
	
    return 1;
}

/***************************************************************************************************
-
-	Euler
-		Do a single Euler Integration Step
-
***************************************************************************************************/
int Mesh::Euler( double* bbox, double* mbb, double endTime, Particle &particle )
{

	//Adjust step to fit in time constraint
	double dt = particle.getStepSize();
    if( particle.t + dt > endTime )
    {
        dt = endTime - particle.t;
		if( dt < epsilon ) return 0;
    }

    double vel[3];
    
    getVelocity( particle, vel );

	double x = particle.x;
	double y = particle.y;
	double z = particle.z;

	//Check if the particle landed in a sink
	if( fabs( vel[0] ) < epsilon && fabs( vel[1] ) < epsilon && fabs( vel[2] ) < epsilon )
	{
		return 0;
	}


	//Move particle
    particle.x += vel[0]*dt;
    particle.y += vel[1]*dt;
    particle.z += vel[2]*dt;    

    particle.t += dt;

	//Check if Steped to end of Time
	if( particle.t >= endTime )
    {
        return 0;
	}

	return checkStep( particle, endTime, bbox, mbb, vel, dt, 0 );
}

int Mesh::RK4( double* bbox, double* mbb, double endTime, Particle &particle )
{

	int toPrint = 0;
/*
	if( almostEqual( particle.x, 1.80909 ) && almostEqual( particle.y, 0.0 ) && almostEqual( particle.z, -0.205637 ) )
	{
		toPrint = 1;
		cerr << "RK:" << endl;
	}
*/
    double k1[3];
    double k2[3];
    double k3[3];
    double k4[3];

    Particle rkpart = particle;
	double dt = particle.getStepSize();

    if( particle.t + dt > endTime )
    {
        dt = endTime - particle.t;
    }

    getVelocity( rkpart, k1 );     

    rkpart.x = particle.x + 0.5*k1[0]*dt;
    rkpart.y = particle.y + 0.5*k1[1]*dt;
    rkpart.z = particle.z + 0.5*k1[2]*dt;

    getVelocity( rkpart, k2 ); 

    rkpart.x = particle.x + 0.5*k2[0]*dt;
    rkpart.y = particle.y + 0.5*k2[1]*dt;
    rkpart.z = particle.z + 0.5*k2[2]*dt;

    getVelocity( rkpart, k3 ); 
    
    rkpart.x = particle.x + k3[0]*dt;
    rkpart.y = particle.y + k3[1]*dt;
    rkpart.z = particle.z + k3[2]*dt;

    getVelocity( rkpart, k4 ); 

    rkpart.x = particle.x;
    rkpart.y = particle.y;
    rkpart.z = particle.z;

	double vel[3] = {	(1.0/6.0)*( k1[0] + 2*k2[0] + 2*k3[0] + k4[0] ), 
						(1.0/6.0)*( k1[1] + 2*k2[1] + 2*k3[1] + k4[1] ),
						(1.0/6.0)*( k1[2] + 2*k2[2] + 2*k3[2] + k4[2] ) };

    particle.x += vel[0]*dt;
    particle.y += vel[1]*dt;
    particle.z += vel[2]*dt;

    particle.t += dt;

	// In Sink
    if( almostEqual( rkpart.x,particle.x) && almostEqual( rkpart.y,particle.y) && almostEqual( rkpart.z,particle.z) )
    {
		return 0;
    }

	//Check if Steped to end of Time
	if( particle.t >= endTime )
    {
        return 0;
	}

	return checkStep( particle, endTime, bbox, mbb, vel, dt, toPrint ); 

}

int Mesh::REV_RK4( double* bbox, double* mbb, double endTime, Particle &particle )
{

	int toPrint = 0;
	if( almostEqual( particle.x, 1.80909 ) && almostEqual( particle.y, 0.0 ) && almostEqual( particle.z, -0.205637 ) )
	{
		toPrint = 1;
	}

    double k1[3];
    double k2[3];
    double k3[3];
    double k4[3];

    Particle rkpart = particle;
	double dt = -1*particle.getStepSize();

    if( particle.t + dt < endTime )
    {
        dt = particle.t - endTime;
    }

    getVelocity( rkpart, k1 );     

    rkpart.x = particle.x + 0.5*k1[0]*dt;
    rkpart.y = particle.y + 0.5*k1[1]*dt;
    rkpart.z = particle.z + 0.5*k1[2]*dt;

    getVelocity( rkpart, k2 ); 

    rkpart.x = particle.x + 0.5*k2[0]*dt;
    rkpart.y = particle.y + 0.5*k2[1]*dt;
    rkpart.z = particle.z + 0.5*k2[2]*dt;

    getVelocity( rkpart, k3 ); 
    
    rkpart.x = particle.x + k3[0]*dt;
    rkpart.y = particle.y + k3[1]*dt;
    rkpart.z = particle.z + k3[2]*dt;

    getVelocity( rkpart, k4 ); 

    rkpart.x = particle.x;
    rkpart.y = particle.y;
    rkpart.z = particle.z;

	double vel[3] = {	(1.0/6.0)*( k1[0] + 2*k2[0] + 2*k3[0] + k4[0] ), 
						(1.0/6.0)*( k1[1] + 2*k2[1] + 2*k3[1] + k4[1] ),
						(1.0/6.0)*( k1[2] + 2*k2[2] + 2*k3[2] + k4[2] ) };

    particle.x += vel[0]*dt;
    particle.y += vel[1]*dt;
    particle.z += vel[2]*dt;

    particle.t += dt;

	// In Sink
    if( almostEqual( rkpart.x,particle.x) && almostEqual( rkpart.y,particle.y) && almostEqual( rkpart.z,particle.z) )
    {
		return 0;
    }

	//Check if Steped to end of Time
	if( particle.t <= endTime )
    {
        return 0;
	}

	return checkStep( particle, endTime, bbox, mbb, vel, dt, toPrint ); 

}

