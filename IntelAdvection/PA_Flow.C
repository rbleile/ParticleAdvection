#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <fstream>

#include "FlowMap.h"
#include "Mesh.h"

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

#define error_tol_m 0.1
#define error_tol 0.00000001
#define step_min  0.000001
#define step_max  0.1

double eps = 0.00001;

class Particle
{

  public:
    double x;
    double y;
    double z;

    double time;
	int face;
    double stepSize;


    Particle() : x(0), y(0), z(0), time(0), face(-1),stepSize(0.01) {};

	void setInfo( FMParticle p )
	{
		x = p.x;
		y = p.y;
		z = p.z;
		time = p.time;
	}

};

/*
 *
 * Idea:
 *        Provide which Dimension will be constant and which will be a range
 *        Define how many particles will be distributed along each range
 *        define x, y, z positions for each particle
 *
*/

// Containes list of particle information
class ParticleList
{
public:
    Particle* particle;

    int numParticles;

	ParticleList( Particle *p, int nump )
	{
		numParticles = nump;
		particle = p;
	}

	ParticleList( FlowMapCell *cell, int face_id )
	{
		numParticles = 4;

		particle = new Particle[ numParticles ];

		for( int i = 0; i < numParticles; i++ )
		{
			particle[i].setInfo( cell->corner_points[ cell->face[face_id].corner_index[i]].in );
		}

	}

    ParticleList( bool x, bool y, bool z, int numP, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax ){

        numParticles = numP;

        particle = new Particle[ numParticles ];

        //What Dimensions will be Constant
        int numVaryingDims = (int)x + (int)y + (int)z;

        if( numVaryingDims == 0 ){ cerr << "Error - Vary Dims Not valid" << endl; exit(0); }

        int numberPerVaryDim;

        if( numVaryingDims == 1)
        {
            numberPerVaryDim = numParticles;
            if( x )
            {
                double step = (xmax-xmin)/numberPerVaryDim;
                for( int i = 0; i < numberPerVaryDim; i++ )
                {
                    particle[i].x = step*i + xmin;
                }
            }
            else
            {
                for( int i = 0; i < numParticles; i++ )
                {
                    particle[i].x = xmin;
                }
            }
    
            if( y )
            {
                double step = (ymax-ymin)/numberPerVaryDim;
                for( int j = 0; j < numberPerVaryDim; j++ ){
                    for( int i = 0; i < numberPerVaryDim; i++ )
                    {
                        particle[i].y = step*i + ymin;
                    }
                }
            }
            else
            {
                for( int i = 0; i < numParticles; i++ )
                {
                    particle[i].y = ymin;
                }
            }
    
            if( z )
            {
                double step = (zmax-zmin)/numberPerVaryDim;
                for( int j = 0; j < numberPerVaryDim; j++ ){
                    for( int i = 0; i < numberPerVaryDim; i++ )
                    {
                        particle[i].z = step*i + zmin;
                    }
                }
            }
            else
            {
                for( int i = 0; i < numParticles; i++ )
                {
                    particle[i].z = zmin;
                }
            }
        }
        else if( numVaryingDims == 2  )
        {
            numberPerVaryDim = (int) sqrt( numParticles );
            if( !x )
            {
                for( int i = 0; i < numParticles; i++ )
                {
                    particle[i].x = xmin;
                }
            }
            else if( !y )
            {
                for( int i = 0; i < numParticles; i++ )
                {
                    particle[i].y = ymin;
                }
            }
            else if( !z )
            {
                for( int i = 0; i < numParticles; i++ )
                {
                    particle[i].z = zmin;
                }
            }

            if( x && y )
            {
                double stepx = (xmax-xmin)/numberPerVaryDim;
                double stepy = (ymax-ymin)/numberPerVaryDim;
				for( int j = 0; j < numberPerVaryDim; j++ ){
	                for( int i = 0; i < numberPerVaryDim; i++ )
		            {
						int id = i + j*numberPerVaryDim;
			            particle[id].x = stepx*i + xmin;
			            particle[id].y = stepy*j + ymin;
				    }
				}
            }
            else if( x & z )
            {
                double stepx = (xmax-xmin)/numberPerVaryDim;
                double stepz = (zmax-zmin)/numberPerVaryDim;
				for( int j = 0; j < numberPerVaryDim; j++ ){
	                for( int i = 0; i < numberPerVaryDim; i++ )
		            {
						int id = i + j*numberPerVaryDim;
			            particle[id].x = stepx*i + xmin;
			            particle[id].z = stepz*j + zmin;
				    }
				}
            }
            else if( y & z )
            {
                double stepy = (ymax-ymin)/numberPerVaryDim;
                double stepz = (zmax-zmin)/numberPerVaryDim;
				for( int j = 0; j < numberPerVaryDim; j++ ){
	                for( int i = 0; i < numberPerVaryDim; i++ )
		            {
						int id = i + j*numberPerVaryDim;
			            particle[id].y = stepy*i + ymin;
			            particle[id].z = stepz*j + zmin;
				    }
				}
            }
        }
        else
        {
            numberPerVaryDim = (int) pow( numParticles, (1.0/3.0));
            double stepx = (xmax-xmin)/numberPerVaryDim;
            double stepy = (ymax-ymin)/numberPerVaryDim;
            double stepz = (zmax-zmin)/numberPerVaryDim;
			for( int k = 0; k < numberPerVaryDim; k++ ){
				for( int j = 0; j < numberPerVaryDim; j++ ){
				  for( int i = 0; i < numberPerVaryDim; i++ )
					{
						int id = i + j*numberPerVaryDim + k*numberPerVaryDim*numberPerVaryDim;
			            particle[id].x = stepx*i + xmin;
			            particle[id].y = stepy*j + ymin;
				        particle[id].z = stepz*k + zmin;
					}
				}
			}
		}
    }

	void SetOUTPoints( FlowMapCell *cell, int face_id )
	{
		
		numParticles = 4;

		for( int i = 0; i < numParticles; i++ )
		{

cerr << "Setting Up: " << particle[i].time << endl;
			cell->corner_points[ cell->face[face_id].corner_index[i]].SetOUTPoint( particle[i].x, particle[i].y, particle[i].z, particle[i].time );
		}

	}
};


///////////////////////////////////////////////////////////////////////////////////////////////////
// Get the instentaneous velocity of a point in the mesh ( Lerp trilinear )
///////////////////////////////////////////////////////////////////////////////////////////////////
void getVelocity ( Particle particle, Mesh* mesh, double *outVelocity )
{

    int pid000[3];

    mesh->getLogicalCellID( particle.x, particle.y, particle.z, pid000 );

	if( pid000[0] == mesh->nx-1 && mesh->isRegular)
	{
	
		if( mesh->dx*(mesh->nx-1) + mesh->x0 == particle.x )
		{
			pid000[0]--;
		}
	}

	if( pid000[1] == mesh->ny-1 && mesh->isRegular )
	{
		if( mesh->dy*(mesh->ny-1) + mesh->y0 == particle.y )
		{
			pid000[1]--;
		}
	}

	if( pid000[2] == mesh->nz-1 && mesh->isRegular)
	{
		if( mesh->dz*(mesh->nz-1) + mesh->z0 == particle.z )
		{
			pid000[2]--;
		}
	}

    if( pid000[0] > mesh->nx-2 || pid000[1] > mesh->ny-2 || pid000[2] > mesh->nz-2 )
    {
		cerr << mesh->nx-2 << " " << mesh->ny-2 << " " << mesh->nz-2 << endl;
        cerr << "Particle: " << particle.x << " " << particle.y << " " << particle.z << endl;
        cerr << "PID:      " << pid000[0] << " " << pid000[1] << " " << pid000[2] << endl;
        cerr << "Why!!" << endl;
        outVelocity[0] = 0;
        outVelocity[1] = 0;
        outVelocity[2] = 0;
        return;
    }

    double x0 = (( mesh->X == NULL) ? mesh->dx*pid000[0] + mesh->x0 : mesh->X[pid000[0]] );
    double y0 = (( mesh->Y == NULL) ? mesh->dy*pid000[1] + mesh->y0 : mesh->Y[pid000[1]] );
    double z0 = (( mesh->Z == NULL) ? mesh->dz*pid000[2] + mesh->z0 : mesh->Z[pid000[2]] );

    double x1 = (( mesh->X == NULL) ? mesh->dx*(pid000[0]+1) + mesh->x0 : mesh->X[pid000[0]+1] );
    double y1 = (( mesh->Y == NULL) ? mesh->dy*(pid000[1]+1) + mesh->y0 : mesh->Y[pid000[1]+1] );
    double z1 = (( mesh->Z == NULL) ? mesh->dz*(pid000[2]+1) + mesh->z0 : mesh->Z[pid000[2]+1] );

    double xd = ( ( particle.x - x0 ) / ( x1 - x0 ) );
    double yd = ( ( particle.y - y0 ) / ( y1 - y0 ) );
    double zd = ( ( particle.z - z0 ) / ( z1 - z0 ) );

    int id000 = pid000[0] + pid000[1]*( mesh->nx  ) + pid000[2]*( mesh->nx )*( mesh->ny );
    int id001 = 3*(id000 + 1);
    int id010 = 3*(id000 + (mesh->nx));
    int id011 = 3*(id000 + (mesh->nx) + 1);

    int id100 = (id000 + (mesh->nx)*(mesh->ny));
    int id101 = 3*(id100 + 1);
    int id110 = 3*(id100 + (mesh->nx));
    int id111 = 3*(id100 + (mesh->nx) + 1);

    id000 *= 3;
    id100 *= 3;

    double c00x = ( ( mesh->V[id000] )   * (1 - xd) ) + ( ( mesh->V[id001])   * xd );
    double c00y = ( ( mesh->V[1+id000] ) * (1 - xd) ) + ( ( mesh->V[1+id001]) * xd );
    double c00z = ( ( mesh->V[2+id000] ) * (1 - xd) ) + ( ( mesh->V[2+id001]) * xd );

    double c10x = ( ( mesh->V[id010] )   * (1 - xd) ) + ( ( mesh->V[id011])   * xd );
    double c10y = ( ( mesh->V[1+id010] ) * (1 - xd) ) + ( ( mesh->V[1+id011]) * xd );
    double c10z = ( ( mesh->V[2+id010] ) * (1 - xd) ) + ( ( mesh->V[2+id011]) * xd );

    double c01x = ( ( mesh->V[id100] )   * (1 - xd) ) + ( ( mesh->V[id101])   * xd );
    double c01y = ( ( mesh->V[1+id100] ) * (1 - xd) ) + ( ( mesh->V[1+id101]) * xd );
    double c01z = ( ( mesh->V[2+id100] ) * (1 - xd) ) + ( ( mesh->V[2+id101]) * xd );

    double c11x = ( ( mesh->V[id110] )   * (1 - xd) ) + ( ( mesh->V[id111])   * xd );
    double c11y = ( ( mesh->V[1+id110] ) * (1 - xd) ) + ( ( mesh->V[1+id111]) * xd );
    double c11z = ( ( mesh->V[2+id110] ) * (1 - xd) ) + ( ( mesh->V[2+id111]) * xd );

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

void stepBackToBoundry( Particle &particle, double* bbox, double *vel, double dt )
{

#if debug
	cerr << endl;

	cerr << "Step Back: " << particle.x << " " << particle.y << " " << particle.z << endl;
	cerr << "Step Too : " << bbox[0] << " " << bbox[1] << " " << bbox[2] << " " <<  bbox[3] << " " << bbox[4] << " " << bbox[5] << endl;
	cerr << "Vel:       " << vel[0] << " " << vel[1] << " " << vel[2] << endl;
	cerr << "dt:        " << dt << endl;
#endif

	double dtx = 0; double dty = 0; double dtz = 0;

	if( particle.x < bbox[0] )
	{
		double xp = bbox[0] - particle.x;
		double vdt = particle.x - ( particle.x - vel[0]*dt);
		dtx = (xp / vdt) * dt;
	}
	else if( particle.x > bbox[1] )
	{
		double xp = bbox[1] - particle.x;
		double vdt = particle.x - ( particle.x - vel[0]*dt);
		dtx = (xp / vdt ) * dt;
	}

	if( particle.y < bbox[2] )
	{
		double yp = bbox[2] - particle.y;
		double vdt = particle.y - ( particle.y - vel[1]*dt);
		dty = (yp / vdt)*dt;
	}
	else if( particle.y > bbox[3] )
	{
		double yp = bbox[3] - particle.y;
		double vdt = particle.y - ( particle.y - vel[1]*dt);
		dty = (yp / vdt)*dt;
	}

	if( particle.z < bbox[4] )
	{
		double zp = bbox[4] - particle.z;
		double vdt = particle.z - ( particle.z - vel[2]*dt);
		dtz = ( zp / vdt ) * dt;
	}
	else if( particle.z > bbox[5] )
	{
		double zp = bbox[5] - particle.z;
		double vdt = particle.z - ( particle.z - vel[2]*dt);
		dtz = ( zp / vdt ) * dt;
	}

	double ndt = 0;

	double adtx = fabs( dtx );
	double adty = fabs( dty );
	double adtz = fabs( dtz );

	if( adtx >= adty && adtx >= adtz )
	{
		ndt = dtx;
	}
	else if( adty >= adtx && adty >= adtz )
	{
		ndt = dty;
	}
	else if( adtz >= adtx && adtz >= adty )
	{
		ndt = dtz;
	}

	particle.x += vel[0]*ndt;
	particle.y += vel[1]*ndt;
	particle.z += vel[2]*ndt;
	particle.time -= fabs( ndt );

#if debug

	cerr << "dts:        " << dtx << " " << dty << " " << dtz << " " << ndt << endl;

	cerr << "atd:        " <<  adtx << " " << adty << " " << adtz << endl;  

	cerr << "Steped Too: " << particle.x << " " << particle.y << " " << particle.z << endl;

cerr << endl;

#endif

}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Check if the particle has left the advection range
///////////////////////////////////////////////////////////////////////////////////////////////////
int checkBoundary( Particle particle, double* bbox )
{

    if( particle.x >= bbox[0] && particle.x <= bbox[1] &&
        particle.y >= bbox[2] && particle.y <= bbox[3] &&
        particle.z >= bbox[4] && particle.z <= bbox[5] )
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Check if the particle has left the advection range
///////////////////////////////////////////////////////////////////////////////////////////////////
int checkOnBoundary( Particle particle, double* bbox )
{

    if( 
		(	( particle.x >= (bbox[0]-eps) && particle.x <= (bbox[0]+eps) )
			|| 
			( particle.x >= (bbox[1]-eps) && particle.x <= (bbox[1]+eps) ) ) 
		||
		(	( particle.y >= (bbox[2]-eps) && particle.y <= (bbox[2]+eps) )
			|| 
			( particle.y >= (bbox[3]-eps) && particle.y <= (bbox[3]+eps) ) ) 
		||
		(	( particle.z >= (bbox[4]-eps) && particle.z <= (bbox[4]+eps) )
			|| 
			( particle.z >= (bbox[5]-eps) && particle.z <= (bbox[5]+eps) ) )
	  )
    {
        return 1;
    }
    else
    {
        return 2;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Do a single advection step in Euler
///////////////////////////////////////////////////////////////////////////////////////////////////
int Euler( Mesh* mesh, Particle &particle, double dt, double endTime, double* bbox )
{
    
    if( particle.time + dt > endTime )
    {
		//cerr << "Step Too large for time: " << particle.time << " " << dt << " " << endTime << endl;
        dt = endTime - particle.time;
    }

    double vel[3];
    
    getVelocity( particle, mesh, vel );

	double pos[3] = { particle.x, particle.y, particle.z };

    particle.x = pos[0] + vel[0]*dt;
    particle.y = pos[1] + vel[1]*dt;
    particle.z = pos[2] + vel[2]*dt;    
    
    particle.time += dt;

    if( particle.time >= endTime )
    {
		cerr << "Step Too large for time: " << particle.time << " " << dt << " " << endTime << endl;
        return 0;
	}else if( vel[0] == 0 && vel[1] == 0 && vel[2] == 0 )
	{
		cerr << "Sink!" << endl;
		return 0;
	}
    else
    {
		int stat = checkBoundary( particle, bbox );
		if( stat ){ 
			stat = checkOnBoundary( particle, bbox ); 
			if( stat == 1 )
			{
				particle.face = -1;
			}
		}
		cerr << "Stat: " << stat << endl;
		if( stat == 0 )
		{
			stepBackToBoundry( particle, bbox, vel, dt );
			particle.face = -1;
			stat = 0;
		}
        return stat;
    }

}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Do a single advection step in RK4
///////////////////////////////////////////////////////////////////////////////////////////////////
int RungeKutta4( Mesh* mesh, Particle &particle, double dt, double endTime )
{

    double k1[3];
    double k2[3];
    double k3[3];
    double k4[3];

    Particle rkpart = particle;

    if( particle.time + dt > endTime )
    {
        dt = endTime - particle.time;
    }

    getVelocity( rkpart, mesh, k1 );     

    rkpart.x = particle.x + 0.5*k1[0]*dt;
    rkpart.y = particle.y + 0.5*k1[1]*dt;
    rkpart.z = particle.z + 0.5*k1[2]*dt;

    getVelocity( rkpart, mesh, k2 ); 

    rkpart.x = particle.x + 0.5*k2[0]*dt;
    rkpart.y = particle.y + 0.5*k2[1]*dt;
    rkpart.z = particle.z + 0.5*k2[2]*dt;

    getVelocity( rkpart, mesh, k3 ); 
    
    rkpart.x = particle.x + k3[0]*dt;
    rkpart.y = particle.y + k3[1]*dt;
    rkpart.z = particle.z + k3[2]*dt;

    getVelocity( rkpart, mesh, k4 ); 

    rkpart.x = particle.x;
    rkpart.y = particle.y;
    rkpart.z = particle.z;

    particle.x += (dt/6.0)*( k1[0] + 2*k2[0] + 2*k3[0] + k4[0] );
    particle.y += (dt/6.0)*( k1[1] + 2*k2[1] + 2*k3[1] + k4[1] );
    particle.z += (dt/6.0)*( k1[2] + 2*k2[2] + 2*k3[2] + k4[2] );

    particle.time += dt;

    if( particle.time >= endTime || rkpart.x == particle.x && rkpart.y == particle.y && rkpart.z == particle.z  )
    {
        return 0;
    }else
    {
        return 1;
    }

}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Do the RKF45 single step calculation
///////////////////////////////////////////////////////////////////////////////////////////////////
double RungeKuttaFehlbergCalculation( Mesh* mesh, Particle &particle, double rk5[3] )
{

    Particle pos = particle;
    double h = particle.stepSize;
    double k1[3];
    double k2[3];
    double k3[3];
    double k4[3];
    double k5[3];
    double k6[3];

//cerr << "Advecting 1 " << pos.x << " " << pos.y << " " << pos.z << endl;
    
    getVelocity( pos, mesh, k1 );
    
    k1[0] *= h;
    k1[1] *= h;
    k1[2] *= h;
    
    pos.x = particle.x + (1.0/4.0) * k1[0];
    pos.y = particle.y + (1.0/4.0) * k1[1];
    pos.z = particle.z + (1.0/4.0) * k1[2];
    
//cerr << "Advecting 2 " << pos.x << " " << pos.y << " " << pos.z <<  endl;
    getVelocity( pos, mesh, k2 );
        
    k2[0] *= h;
    k2[1] *= h;
    k2[2] *= h;
    
    pos.x = particle.x + (3.0/32.0)*k1[0] + (9.0/32.0)*k2[0];
    pos.y = particle.y + (3.0/32.0)*k1[1] + (9.0/32.0)*k2[1];
    pos.z = particle.z + (3.0/32.0)*k1[2] + (9.0/32.0)*k2[2];
    
//cerr << "Advecting 3 " << pos.x << " " << pos.y << " " << pos.z <<  endl;
    getVelocity( pos, mesh, k3 );
    
    k3[0] *= h;
    k3[1] *= h;
    k3[2] *= h;
    
    pos.x = particle.x + (1932.0/2197.0)*k1[0] - ( 7200.0 / 2197.0 )*k2[0] + (7296.0/2197.0)*k3[0];
    pos.y = particle.y + (1932.0/2197.0)*k1[1] - ( 7200.0 / 2197.0 )*k2[1] + (7296.0/2197.0)*k3[1];
    pos.z = particle.z + (1932.0/2197.0)*k1[2] - ( 7200.0 / 2197.0 )*k2[2] + (7296.0/2197.0)*k3[2];
        
//cerr << "Advecting 4 " << pos.x << " " << pos.y << " " << pos.z << endl;
    getVelocity( pos, mesh, k4 );
    
    k4[0] *= h;
    k4[1] *= h;
    k4[2] *= h;
    
    pos.x = particle.x + (439.0/216.0)*k1[0] - 8.0*k2[0] + (3680.0/513.0)*k3[0] - (845.0/4104.0)*k4[0];
    pos.y = particle.y + (439.0/216.0)*k1[1] - 8.0*k2[1] + (3680.0/513.0)*k3[1] - (845.0/4104.0)*k4[1];
    pos.z = particle.z + (439.0/216.0)*k1[2] - 8.0*k2[2] + (3680.0/513.0)*k3[2] - (845.0/4104.0)*k4[2];
    
//cerr << "Advecting 5 " << pos.x << " " << pos.y << " " << pos.z << endl;
    getVelocity( pos, mesh, k5 );
    
    k5[0] *= h;
    k5[1] *= h;
    k5[2] *= h;
        
    pos.x = particle.x - (8.0/27.0)*k1[0] + 2.0*k2[0] - (3544.0/2565.0)*k3[0] + (1859.0/4104.0)*k4[0] - (11.0/40.0)*k5[0];
    pos.y = particle.y - (8.0/27.0)*k1[1] + 2.0*k2[1] - (3544.0/2565.0)*k3[1] + (1859.0/4104.0)*k4[1] - (11.0/40.0)*k5[1];
    pos.z = particle.z - (8.0/27.0)*k1[2] + 2.0*k2[2] - (3544.0/2565.0)*k3[2] + (1859.0/4104.0)*k4[2] - (11.0/40.0)*k5[2];
    
//cerr << "Advecting 6 " << pos.x << " " << pos.y << " " << pos.z <<  endl;
    getVelocity( pos, mesh, k6 );
    
    k6[0] *= h;
    k6[1] *= h;
    k6[2] *= h;

    double rk4[3]; 

    rk4[0] = particle.x;
    rk5[0] = particle.x;
    rk4[1] = particle.y;
    rk5[1] = particle.y;
    rk4[2] = particle.z;
    rk5[2] = particle.z;

    for( int i = 0; i < 3; i++)
    { 
        rk4[i] += (25.0/216.0)*k1[i] + (1408.0/2565.0)*k3[i] + (2197.0/4104.0)*k4[i] - (1.0/5.0)*k5[i];
        rk5[i] += (16.0/135.0)*k1[i] + (6656.0/12825.0)*k3[i] + (28561.0/56430.0)*k4[i] - (9.0/50.0)*k5[i] + (2.0/55.0)*k6[i];    
    }

//cerr << "Advecting 7 " << rk4[0] << " " << rk4[1] << " " << rk4[2] <<  " " << rk5[0] << " " << rk5[1] << " " << rk5[2] <<  endl;

    return sqrt( pow( rk5[0] - rk4[0], 2) +  pow( rk5[1] - rk4[1], 2) +  pow( rk5[2] - rk4[2], 2) );
   

}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Do a single RKF45 step -> advect or change step size 
///////////////////////////////////////////////////////////////////////////////////////////////////
int RungeKuttaFehlbergStep( Mesh* mesh, Particle &particle, double endTime )
{

    double timeLeft = endTime - particle.time;

    if( timeLeft == 0 )
    {
        return 0;
    }else if ( particle.stepSize > timeLeft ) 
    {
        particle.stepSize = timeLeft;
    }

    double tempPos[3];

    double error_45 = RungeKuttaFehlbergCalculation( mesh, particle, tempPos );

//cerr << error_45 << endl;

    if( error_tol * error_tol_m <= error_45 && error_45 <= error_tol )
    {
        particle.x = tempPos[0];
        particle.y = tempPos[1];
        particle.z = tempPos[2];
        particle.time += particle.stepSize;
        return 1;
    } 
    else
    {
        if( error_45 < error_tol * error_tol_m && particle.stepSize*2.0 <= step_max )
        {
            particle.x = tempPos[0];
            particle.y = tempPos[1];
            particle.z = tempPos[2];
            particle.time += particle.stepSize;
            
            particle.stepSize *= 2.0; //Step Size Too Small but don't throw away calculation
            //particle.stepSize *= (pow( ( (error_tol * particle.stepSize) / (2*error_45) ) , (1.0/4.0)));

            return 1;
        }
        else if( error_45 > error_tol && particle.stepSize/2 >= step_min )
        {
            particle.stepSize *= (pow( ( (error_tol * particle.stepSize) / (2*error_45) ) , (1.0/4.0)));
            //particle.stepSize /= 2.0;
/*
            double s = 0.840896 * pow(( (error_max * particle.stepSize) / error_45 ), (1.0/4.0));
            particle.stepSize *= s;
*/
            //cerr << "Change Step Size: 2 /" << endl;
            return 2;
        }
        else
        {
        // Forced to keep the current step size regardless of error and advect particle
        //     char* errorString[128]; 
        //    sprintf( errorString, "Error Exceeds Limits (%e) and StepSize at bounds (%e)", error_45, stepSize );
            
            particle.x = tempPos[0];
            particle.y = tempPos[1];
            particle.z = tempPos[2];
            particle.time += particle.stepSize;
			//cerr << "Change Position: Error" << endl;
			return 1;
		}

    }

}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Do a single advection step for one particle for a given method
///////////////////////////////////////////////////////////////////////////////////////////////////
int advectParticle( Mesh* mesh, Particle &particle, double *bbox, double endTime )
{

    int status = 1;

    #ifdef RK4
 
        /* Do Runge Kutta 4  */
        double dt = particle.stepSize;
        status = RungeKutta4( mesh, particle, dt, endTime );
        /* status: 0 = terminated / 1 = advected */


    #elif RKF45

        /* Do Runge Kutta Fehlberg 4-5 */

        status = RungeKuttaFehlbergStep( mesh, particle, endTime );
        /* status: 0 = terminated / 1 = advected / 2 = stepSize changed  */

    #else

        /* Do Euler */

        double dt = particle.stepSize;
        status = Euler( mesh, particle, dt, endTime, bbox );

		//cerr << "Status: " << status << endl;

    #endif

        return status;
}

/*Track Particles using Eulerian Method*/
void trackParticles( Mesh* mesh, ParticleList* particleList, double* bbox, double endTime )
{

    int numP = particleList->numParticles;

    for( int i = 0; i < numP; i++ )
    {

        int status = 1;

        Particle &particle = particleList->particle[i];

		cerr << i << " BB: " << bbox[0] << " " << bbox[1] << " " << bbox[2] << " " <<  bbox[3] << " " << bbox[4] << " " << bbox[5] << endl;

        while( status )
        {

			cerr << i << ": " << particle.x << " " << particle.y << " " << particle.z << endl; 

            status = advectParticle( mesh, particle, bbox, endTime );

		cerr << "particle: " << i << "\t" << status << endl;

            if( status )
            {
                status = checkBoundary( particle, bbox );
            }

        }
    }
}

int computeCellFace( Particle &particle, Mesh* mesh, double* bbox, double* out_bbox, int &cellID )
{

//See if particle is in the mesh
	int test = checkBoundary( particle, bbox );
	if( test == 0 )
	{
		
		return 0;
	}

	//Get a second particle to advect a small step and determing the output cell that should used in the event of no acceptable faces are found
	Particle outParticle = particle;

	Euler( mesh, outParticle, outParticle.stepSize, 1, bbox );

//If so find all cells that contain the particle
	for( int i = 0; i < mesh->num_cells; i++ )
	{
		FlowMapCell* fmc = &mesh->fmc[i];
		
		test = checkBoundary( outParticle, fmc->bbox );
		if( test )
		{
			//Set this just in case the faces have not accepted flows
			cellID = i;
			out_bbox = fmc->bbox;
		}

		test = checkBoundary( particle, fmc->bbox );
		if( test )
		{
//Then decide if particle is on a face
			test = checkOnBoundary( particle, fmc->bbox );
			if( test == 1 )
			{
			
//If so gather which face or set of faces
				int face_id;
				int edge_id;
				int corn_id;
				fmc->getCellFace( particle.x, particle.y, particle.z, face_id, edge_id, corn_id );

				if( face_id != -1 )
				{
//return face if particle on a face and accepted flow if not try a new cell
					cerr << "Particle on Face: " << face_id << endl;
					particle.face = face_id;
					if( fmc->face[ face_id].accepted_flow )
					{
						cerr << "Face is accepted flow - FLOW" << endl;
						cellID = i;
						return 1;
					}
					else
					{
						cerr << "Face is NOT accepted flow - EULER" << endl;
						cerr << "Try another cell" << endl;
						continue;
					}
				}
				else
				{
//if face not get a list of faces from edges or corners
					int* adjacentFaces;
					int numAdjacentFaces;

					if( edge_id != -1  )
					{
						numAdjacentFaces=2;
						adjacentFaces = new int [numAdjacentFaces];
						
						for( int j = 0; j < 2; j++ )
						{
							adjacentFaces[j] = EdgesToFaces[edge_id][j];
						}
					}
					else
					{
						if( corn_id != -1 )
						{
							numAdjacentFaces=2;
							adjacentFaces = new int [numAdjacentFaces];
							
							for( int j = 0; j < 2; j++ )
							{
								adjacentFaces[j] = EdgesToFaces[edge_id][j];
							}
							
						}
						else
						{
							//no matching face or edge or corner
							//Should Not reach this point
							cerr << "Should Not Reach This Condition - ERROR" << endl;
							return 0;
						}
					}

//with list of faces check if any of these faces has an accepted flow. If so that is the face to use for flow
//if not then try a new cell
					for( int j = 0; j < numAdjacentFaces; j++ )
					{
						if( fmc->face[ adjacentFaces[j] ].accepted_flow )
						{
							particle.face = adjacentFaces[j];
							cellID = i;
							return 1;
						}
					}

				}

			}
			else
			{
				//Particle is internal to a cell and theirfore not on a face
				cellID = i;
				out_bbox = fmc->bbox;
				return 2;
			}
		}
	}

	//Particle face with accepted flow not found
	
	particle.face = -1;
	return 2;
}

/* Track Particles using FLow Maps when possible */
void trackParticlesFlow( Mesh *mesh, ParticleList* particleList, double* bbox, double endTime )
{
	int numP = particleList->numParticles;

	int numEuler = 0;
	int numFlow = 0;

	for( int i = 0; i < numP; i++ )
	{
		int status = 1;
		Particle &particle = particleList->particle[i];

		double cellBB[6];
		int cellID;

		while( status )
		{
			cerr << "Particle ID: " << i << "Time: " << particle.time << "\tX: " << particle.x << "\tY: " << particle.y << "\tZ: " << particle.z << "\tFace: " << particle.face << "\tStatus: " << status << "\tNumEuler: " << numEuler << "\tNumFlow " << numFlow << endl;

			if( particle.face == -1)
			{
				status = computeCellFace( particle, mesh, bbox, cellBB, cellID );
			}

			if( status == 0 )
			{
				break;
			}

			/*Euler if 2, Try Flow if 1 */ 
			if( status == 2 )
			{
				status = advectParticle( mesh, particle, cellBB, endTime );
				numEuler++;
			}
			else
			{
				int cell_change = 0;
				status = mesh->fmc[cellID].advectOnFlow( particle.x, particle.y, particle.z, particle.time, particle.face, cell_change  );

				if( status == 0)
				{
					cerr << "FLOW: No time change in flow advection - ERROR ?!?" << endl;
				}
				else if( status == 1 )
				{
					cerr << "FLOW: successful advection on flow" << endl;
				}
				else if( status == 2 )
				{
					particle.face = -1;
					cerr << "FLOW: cannot advect must euler instead" << endl;
				}

				if( status != 2 )
					numFlow++;

				if( particle.time > endTime )
				{
					status = 0;
				}
			}
		}
	}

}

void outputParticles( ParticleList *particlesList, char* name)
{

    char buffer[256];
    
#ifdef RK4
    const char algorithm[] = "RK4";
#elif RKF45
    const char algorithm[] = "RKF45";
#else
    const char algorithm[] = "EULER";
#endif

    const char end[5] = ".dat";
    sprintf(buffer, "output%s/Particles_%s%s", algorithm,name,end);

    ofstream particle_output ( buffer, std::ios::out );
    particle_output.precision( 8 );

    for( int i = 0; i < particlesList->numParticles; i++ )
    {
        particle_output << particlesList->particle[i].x << "\t" << particlesList->particle[i].y << "\t" << particlesList->particle[i].z  << "\t" << particlesList->particle[i].time<< endl; 
    }

    particle_output.close();

}


void buildRectilinearMesh()
{

}

void buildRegularMesh()
{

}


int main()
{

    cerr << " START " << endl;


/*
 *
 * Define the Mesh we will be doing advection on
 *
*/

// Number of Points in the XYZ dimensions
    const int nx = 11;
    const int ny = 11;
    const int nz = 11;

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

    cerr << "make in vels" << endl;

    for( int z = 0; z < nz; z++ ){
        for( int y = 0; y < ny; y++ ){
            for( int x = 0; x < nx; x++ )
            {


                double xp = xmin + x*dx;
                double yp = ymin + y*dy;
                double zp = zmin + z*dz;

                double r = sqrt( xp*xp + yp*yp );
                double t = atan( yp/(xp+0.01));
                double px = r * cos( t );
                double py = r * sin( t );
				t = atan( zp / ( yp+0.01) );
				double pz = 3 * r * sin( t );

                int id = 3*(x + y*nx + z*nx*ny);
                
                if( px == 0 && py == 0 )
                {
                    px = 1;
                    py = 1;
                }

                v_field[  id  ] = px;
                v_field[ id+1 ] = py;
                v_field[ id+2 ] = pz;


/*
                int id = 3*(x + y*nx + z*nx*ny);
                v_field[  id  ] = 1.0;
                v_field[ id+1 ] = 0.0;
                v_field[ id+2 ] = 0.0;
*/				
            }
        }
    }


    bool isReg = true;

    Mesh mesh( isReg, nx, ny, nz, v_field, dx, dy, dz, xmin, ymin, zmin );

/*
 *
 * Define the Flow Map we will be laying over the mesh
 *
*/

	int num_cells = (nx-1)*(ny-1)*(nz-1);

	mesh.num_cells = num_cells;
	mesh.fmc = new FlowMapCell [ num_cells ];

	fillCellParticlesIN( &mesh );

	PrintFlowMap( &mesh, cout );

/* Advect Entire mesh of particles to determine the usefulness of the flow map */

	int num_computable = 0;
	int total_faces = 0;

	for( int cell = 0; cell < num_cells; cell++ )
	{
		cerr << "Cell: " << cell << endl;
	
		FlowMapCell &FMcell = mesh.fmc[cell];

		double *bb = FMcell.bbox;

		for( int f = 0; f < 6; f++ )
		{
			cerr << "Face: " << f << endl;

			ParticleList pl( &FMcell, f );

			trackParticles( &mesh, &pl, bb, 1000); 

			cerr << "Out: " << endl;

			pl.SetOUTPoints( &FMcell, f );	
			
		}	

		FMcell.computeOUTCellFaces();

		total_faces += 6;
		num_computable += FMcell.computeAcceptableFaces();

	}

	cerr << "Number of Computable Faces: " << num_computable << endl;
	cerr << "Number of total  Faces:     " << total_faces << endl;

	cerr << "Percentage of computable faces: " << ((double)num_computable / (double)total_faces ) * 100 << endl;

	
	PrintFlowMap( &mesh, cout );
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
}
