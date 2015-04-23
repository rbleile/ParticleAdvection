#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <fstream>

#include "FlowMap.h"
#include "Mesh.h"

using std::cerr;
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

class Particle
{

  public:
    float x;
    float y;
    float z;

    float time;
    float stepSize;


    Particle() : x(0), y(0), z(0), time(0), stepSize(0.01) {};

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

	ParticleList( FlowMapCell *cell, int face_id )
	{
		numParticles = 4;

		particle = new Particle[ numParticles ];

		for( int i = 0; i < numParticles; i++ )
		{
			particle[i].setInfo( cell->corner_points[ cell->face[face_id].corner_index[i]].in );
		}

	}

    ParticleList( bool x, bool y, bool z, int numP, float xmin, float xmax, float ymin, float ymax, float zmin, float zmax ){

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
                float step = (xmax-xmin)/numberPerVaryDim;
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
                float step = (ymax-ymin)/numberPerVaryDim;
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
                float step = (zmax-zmin)/numberPerVaryDim;
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
                float stepx = (xmax-xmin)/numberPerVaryDim;
                float stepy = (ymax-ymin)/numberPerVaryDim;
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
                float stepx = (xmax-xmin)/numberPerVaryDim;
                float stepz = (zmax-zmin)/numberPerVaryDim;
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
                float stepy = (ymax-ymin)/numberPerVaryDim;
                float stepz = (zmax-zmin)/numberPerVaryDim;
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
            float stepx = (xmax-xmin)/numberPerVaryDim;
            float stepy = (ymax-ymin)/numberPerVaryDim;
            float stepz = (zmax-zmin)/numberPerVaryDim;
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

			cell->corner_points[ cell->face[face_id].corner_index[i]].SetOUTPoint( particle[i].x, particle[i].y, particle[i].z, particle[i].time );
		}

	}
};


///////////////////////////////////////////////////////////////////////////////////////////////////
// Get the instentaneous velocity of a point in the mesh ( Lerp trilinear )
///////////////////////////////////////////////////////////////////////////////////////////////////
void getVelocity ( Particle particle, Mesh* mesh, float *outVelocity )
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

    float x0 = (( mesh->X == NULL) ? mesh->dx*pid000[0] + mesh->x0 : mesh->X[pid000[0]] );
    float y0 = (( mesh->Y == NULL) ? mesh->dy*pid000[1] + mesh->y0 : mesh->Y[pid000[1]] );
    float z0 = (( mesh->Z == NULL) ? mesh->dz*pid000[2] + mesh->z0 : mesh->Z[pid000[2]] );

    float x1 = (( mesh->X == NULL) ? mesh->dx*(pid000[0]+1) + mesh->x0 : mesh->X[pid000[0]+1] );
    float y1 = (( mesh->Y == NULL) ? mesh->dy*(pid000[1]+1) + mesh->y0 : mesh->Y[pid000[1]+1] );
    float z1 = (( mesh->Z == NULL) ? mesh->dz*(pid000[2]+1) + mesh->z0 : mesh->Z[pid000[2]+1] );

    float xd = ( ( particle.x - x0 ) / ( x1 - x0 ) );
    float yd = ( ( particle.y - y0 ) / ( y1 - y0 ) );
    float zd = ( ( particle.z - z0 ) / ( z1 - z0 ) );

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

    float c00x = ( ( mesh->V[id000] )   * (1 - xd) ) + ( ( mesh->V[id001])   * xd );
    float c00y = ( ( mesh->V[1+id000] ) * (1 - xd) ) + ( ( mesh->V[1+id001]) * xd );
    float c00z = ( ( mesh->V[2+id000] ) * (1 - xd) ) + ( ( mesh->V[2+id001]) * xd );

    float c10x = ( ( mesh->V[id010] )   * (1 - xd) ) + ( ( mesh->V[id011])   * xd );
    float c10y = ( ( mesh->V[1+id010] ) * (1 - xd) ) + ( ( mesh->V[1+id011]) * xd );
    float c10z = ( ( mesh->V[2+id010] ) * (1 - xd) ) + ( ( mesh->V[2+id011]) * xd );

    float c01x = ( ( mesh->V[id100] )   * (1 - xd) ) + ( ( mesh->V[id101])   * xd );
    float c01y = ( ( mesh->V[1+id100] ) * (1 - xd) ) + ( ( mesh->V[1+id101]) * xd );
    float c01z = ( ( mesh->V[2+id100] ) * (1 - xd) ) + ( ( mesh->V[2+id101]) * xd );

    float c11x = ( ( mesh->V[id110] )   * (1 - xd) ) + ( ( mesh->V[id111])   * xd );
    float c11y = ( ( mesh->V[1+id110] ) * (1 - xd) ) + ( ( mesh->V[1+id111]) * xd );
    float c11z = ( ( mesh->V[2+id110] ) * (1 - xd) ) + ( ( mesh->V[2+id111]) * xd );

    float c0x = c00x * ( 1 - yd ) + c10x * yd;
    float c0y = c00y * ( 1 - yd ) + c10y * yd;
    float c0z = c00z * ( 1 - yd ) + c10z * yd;

    float c1x = c01x * ( 1 - yd ) + c11x * yd;
    float c1y = c01y * ( 1 - yd ) + c11y * yd;
    float c1z = c01z * ( 1 - yd ) + c11z * yd;

    outVelocity[0] = c0x * ( 1 - zd ) + c1x * zd;    
    outVelocity[1] = c0y * ( 1 - zd ) + c1y * zd;    
    outVelocity[2] = c0z * ( 1 - zd ) + c1z * zd;    

}

void stepBackToBoundry( Particle &particle, float* bbox, float *vel, float dt )
{

#if debug
	cerr << endl;

	cerr << "Step Back: " << particle.x << " " << particle.y << " " << particle.z << endl;
	cerr << "Step Too : " << bbox[0] << " " << bbox[1] << " " << bbox[2] << " " <<  bbox[3] << " " << bbox[4] << " " << bbox[5] << endl;
	cerr << "Vel:       " << vel[0] << " " << vel[1] << " " << vel[2] << endl;
	cerr << "dt:        " << dt << endl;
#endif

	float dtx = 0; float dty = 0; float dtz = 0;

	if( particle.x < bbox[0] )
	{
		float xp = bbox[0] - particle.x;
		float vdt = particle.x - ( particle.x - vel[0]*dt);
		dtx = (xp / vdt) * dt;
	}
	else if( particle.x > bbox[1] )
	{
		float xp = bbox[1] - particle.x;
		float vdt = particle.x - ( particle.x - vel[0]*dt);
		dtx = (xp / vdt ) * dt;
	}

	if( particle.y < bbox[2] )
	{
		float yp = bbox[2] - particle.y;
		float vdt = particle.y - ( particle.y - vel[1]*dt);
		dty = (yp / vdt)*dt;
	}
	else if( particle.y > bbox[3] )
	{
		float yp = bbox[3] - particle.y;
		float vdt = particle.y - ( particle.y - vel[1]*dt);
		dty = (yp / vdt)*dt;
	}

	if( particle.z < bbox[4] )
	{
		float zp = bbox[4] - particle.z;
		float vdt = particle.z - ( particle.z - vel[2]*dt);
		dtz = ( zp / vdt ) * dt;
	}
	else if( particle.z > bbox[5] )
	{
		float zp = bbox[5] - particle.z;
		float vdt = particle.z - ( particle.z - vel[2]*dt);
		dtz = ( zp / vdt ) * dt;
	}

	float ndt = 0;

	float adtx = fabs( dtx );
	float adty = fabs( dty );
	float adtz = fabs( dtz );

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
int checkBoundary( Particle particle, float* bbox )
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
// Do a single advection step in Euler
///////////////////////////////////////////////////////////////////////////////////////////////////
int Euler( Mesh* mesh, Particle &particle, float dt, float endTime, float* bbox )
{
    
    if( particle.time + dt > endTime )
    {
		//cerr << "Step Too large for time: " << particle.time << " " << dt << " " << endTime << endl;
        dt = endTime - particle.time;
    }

    float vel[3];
    
    getVelocity( particle, mesh, vel );

	float pos[3] = { particle.x, particle.y, particle.z };

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
		//cerr << "Stat: " << stat << endl;
		if( stat == 0 )
		{
			stepBackToBoundry( particle, bbox, vel, dt );
		}
        return stat;
    }

}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Do a single advection step in RK4
///////////////////////////////////////////////////////////////////////////////////////////////////
int RungeKutta4( Mesh* mesh, Particle &particle, float dt, float endTime )
{

    float k1[3];
    float k2[3];
    float k3[3];
    float k4[3];

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
float RungeKuttaFehlbergCalculation( Mesh* mesh, Particle &particle, float rk5[3] )
{

    Particle pos = particle;
    float h = particle.stepSize;
    float k1[3];
    float k2[3];
    float k3[3];
    float k4[3];
    float k5[3];
    float k6[3];

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

    float rk4[3]; 

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
int RungeKuttaFehlbergStep( Mesh* mesh, Particle &particle, float endTime )
{

    float timeLeft = endTime - particle.time;

    if( timeLeft == 0 )
    {
        return 0;
    }else if ( particle.stepSize > timeLeft ) 
    {
        particle.stepSize = timeLeft;
    }

    float tempPos[3];

    float error_45 = RungeKuttaFehlbergCalculation( mesh, particle, tempPos );

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
            float s = 0.840896 * pow(( (error_max * particle.stepSize) / error_45 ), (1.0/4.0));
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
int advectParticle( Mesh* mesh, Particle &particle, float *bbox, float endTime )
{

    int status = 1;

    #ifdef RK4
 
        /* Do Runge Kutta 4  */
        float dt = particle.stepSize;
        status = RungeKutta4( mesh, particle, dt, endTime );
        /* status: 0 = terminated / 1 = advected */


    #elif RKF45

        /* Do Runge Kutta Fehlberg 4-5 */

        status = RungeKuttaFehlbergStep( mesh, particle, endTime );
        /* status: 0 = terminated / 1 = advected / 2 = stepSize changed  */

    #else

        /* Do Euler */

        float dt = particle.stepSize;
        status = Euler( mesh, particle, dt, endTime, bbox );

		//cerr << "Status: " << status << endl;

    #endif

        return status;
}



/*
 *
 * Stephanie - here
 *
 * */

void trackParticles( Mesh* mesh, ParticleList* particleList, float* bbox, float endTime )
{

    int numP = particleList->numParticles;

    for( int i = 0; i < numP; i++ )
    {

        int status = 1;

        Particle &particle = particleList->particle[i];

		//cerr << i << " BB: " << bbox[0] << " " << bbox[1] << " " << bbox[2] << " " <<  bbox[3] << " " << bbox[4] << " " << bbox[5] << endl;

        while( status )
        {

			//cerr << i << ": " << particle.x << " " << particle.y << " " << particle.z << endl; 

            status = advectParticle( mesh, particle, bbox, endTime );

			//cerr << "particle: " << i << "\t" << status << endl;

/*
            if( status == 1 )
            {
                status = checkBoundary( particle, bbox );

				if( status == 0 )
				{
					
				}
            }
*/
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


void advectParticleOnFace()
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

    float* v_field = new float [ 3 * nx * ny * nz ];

// The minimum Corner of the Mesh
    float xmin = -10.0;
    float ymin = -10.0;
    float zmin = -10.0;

// The Maximum Corner of the Mesh
    float xmax = 10.0;
    float ymax = 10.0;
    float zmax = 10.0;

// Step Size on Mesh
    float dx = (xmax-xmin)/(nx-1);
    float dy = (ymax-ymin)/(ny-1);
    float dz = (zmax-zmin)/(nz-1);

    cerr << "makein vels" << endl;

    for( int z = 0; z < nz; z++ ){
        for( int y = 0; y < ny; y++ ){
            for( int x = 0; x < nx; x++ )
            {


                float xp = xmin + x*dx;
                float yp = ymin + y*dy;
                float zp = zmin + z*dz;

                float r = sqrt( xp*xp + yp*yp );
                float t = atan( yp/(xp+0.01));
                float px = r * cos( t );
                float py = r * sin( t );
				t = atan( zp / ( yp+0.01) );
				float pz = 3 * r * sin( t );

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

	PrintFlowMap( &mesh, cerr );

/* Advect Entire mesh of particles to determine the usefulness of the flow map */

	int num_computable = 0;
	int total_faces = 0;

	for( int cell = 0; cell < num_cells; cell++ )
	{
	
		FlowMapCell &FMcell = mesh.fmc[cell];

		float *bb = FMcell.bbox;

		for( int f = 0; f < 6; f++ )
		{

			ParticleList pl( &FMcell, f );

			trackParticles( &mesh, &pl, bb, 1000); 

			pl.SetOUTPoints( &FMcell, f );	
			
		}	

		FMcell.computeOUTCellFaces();

		total_faces += 6;
		num_computable += FMcell.computeAcceptableFaces();

	}
	PrintFlowMap( &mesh, cerr );

	cerr << "Number of Computable Faces: " << num_computable << endl;
	cerr << "Number of total  Faces:     " << total_faces << endl;

	cerr << "Percentage of computable faces: " << ((float)num_computable / (float)total_faces ) * 100 << endl;


}
