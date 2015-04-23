#ifndef FLOWMAP_H
#define FLOWMAP_H

#include <ostream>

using std::ostream;
using std::cerr;
using std::endl;

const int PointsToFaces[8][3] = 
	{
		{0,3,4},
		{0,1,4},
		{2,3,4},
		{1,2,4},
		{0,3,5},
		{0,1,5},
		{2,3,5},
		{1,2,5}
	};

const int EdgesToFaces[12][2] = 
	{
		{0, 4},
		{1, 4},
		{2, 4},
		{3, 4},
		{0, 5},
		{1, 5},
		{2, 5},
		{3, 5},
		{0, 3},
		{0, 1},
		{1, 2},
		{2, 3},
	};

const int FaceToFace[6]= { 2, 3, 0, 1, 5, 4 };

/**************************************************************************************************
 *
 *	FMParticle
 *		Flow Map Particle Information. Associated with the advected corners of a face to hold 
 *		input and output informaiton for out lagrangian flow particles
 *
 *
**************************************************************************************************/
class FMParticle
{
  public:
	double x;
	double y;
	double z;
	double time;
	int face_index;
	int corner_index;
	int edge_index;

	FMParticle();
	void printParticle( ostream &stream );


};

FMParticle::FMParticle()
{

	x = 0;
	y = 0;
	z = 0;
	time = 0;
	face_index = -1;
	corner_index = -1;
	edge_index = -1;

}

double epsilon = 0.00001;

inline bool operator==( const FMParticle &p1, const FMParticle &p2 )
{

	if(
		( p1.x <= (p2.x+epsilon) && p1.x >= (p2.x-epsilon) )  &&
		( p1.y <= (p2.y+epsilon) && p1.y >= (p2.y-epsilon) )  &&
		( p1.z <= (p2.z+epsilon) && p1.z >= (p2.z-epsilon) )
	  )
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool operator!=( const FMParticle &p1, const FMParticle &p2 )
{
	return !(p1==p2);
}

void FMParticle::printParticle( ostream &stream )
{
	stream << "\t\t\t< X: " << x << ", Y: " << y << ", Z: " << z << ", T: " << time << ", FI: " << face_index << ", CI: " << corner_index << ", EI: " << edge_index  << " >" << endl;
}

/**************************************************************************************************
 *
 *	FMCornerPoint
 *		Flow Map Corners of the faces. Holds the particle informaiton associated with each corner
 *		of each face. IN is particle position and time at the corner. Out is the advected particle 
 *		after having been advected from the in location in the flow map cell
 *
 *
**************************************************************************************************/
class FMCornerPoint
{
  public:
	FMParticle in;
	FMParticle out;

	void SetINPoint( double xp, double yp, double zp );
    void SetINCornerID( int id );
	void SetOUTPoint( double xp, double yp, double zp );
	void SetOUTPoint( double xp, double yp, double zp, double time );
    void SetOUTCornerID( int id );
    void SetOUTCornerFace( int id );

	void printCorner( ostream &stream );

};

void FMCornerPoint::SetINPoint( double xp, double yp, double zp )
{
	in.x = xp;
	in.y = yp;
	in.z = zp;
}
void FMCornerPoint::SetINCornerID( int id )
{
	in.corner_index = id;
}

void FMCornerPoint::SetOUTPoint( double xp, double yp, double zp )
{
	out.x = xp;	
	out.y = yp;	
	out.z = zp;	
}

void FMCornerPoint::SetOUTPoint( double xp, double yp, double zp, double time )
{
	out.x = xp;	
	out.y = yp;	
	out.z = zp;	
	out.time = time;
}

void FMCornerPoint::SetOUTCornerID( int id )
{
	out.corner_index = id;
}

void FMCornerPoint::SetOUTCornerFace( int id )
{
	out.face_index = id;
}


void FMCornerPoint::printCorner( ostream &stream )
{
	stream << "\t\t" << "Particle IN: " << endl;
	in.printParticle( stream );
	stream << "\t\t" << "Particle OUT: " << endl;
	out.printParticle( stream );
}

/**************************************************************************************************
 *
 *	FMFace
 *		Flow Map cell faces. Each face contains four corner points and determines wherether or not 
 *		advection has been done yet or can be done on this cell
 *
 *
**************************************************************************************************/
class FMFace
{
  public:
	bool accepted_flow;
	int index;
	int corner_index[4];
	int edge_index[4];

	void INIT( int id );
	void SetCorners( );
	void printFace( ostream &stream );

};
void FMFace::INIT( int id )
{
	accepted_flow = false;
	index = id;
	SetCorners();
}

void FMFace::SetCorners()
{
	if( index == 0 )
	{
		corner_index[0] = 0;
		corner_index[1] = 1;
		corner_index[2] = 4;
		corner_index[3] = 5;
		edge_index[0] = 0;
		edge_index[1] = 4;
		edge_index[2] = 8;
		edge_index[3] = 9;
	}
	if( index == 1 )
	{
		corner_index[0] = 1;
		corner_index[1] = 3;
		corner_index[2] = 5;
		corner_index[3] = 7;
		edge_index[0] = 1;
		edge_index[1] = 5;
		edge_index[2] = 9;
		edge_index[3] = 10;
	}
	if( index == 2 )
	{
		corner_index[0] = 2;
		corner_index[1] = 3;
		corner_index[2] = 6;
		corner_index[3] = 7;
		edge_index[0] = 2;
		edge_index[1] = 6;
		edge_index[2] = 10;
		edge_index[3] = 11;
	}
	if( index == 3 )
	{
		corner_index[0] = 0;
		corner_index[1] = 2;
		corner_index[2] = 4;
		corner_index[3] = 6;
		edge_index[0] = 3;
		edge_index[1] = 7;
		edge_index[2] = 11;
		edge_index[3] = 8;
	}
	if( index == 4 )
	{
		corner_index[0] = 0;
		corner_index[1] = 1;
		corner_index[2] = 2;
		corner_index[3] = 3;
		edge_index[0] = 0;
		edge_index[1] = 1;
		edge_index[2] = 2;
		edge_index[3] = 3;
	}
	if( index == 5 )
	{
		corner_index[0] = 4;
		corner_index[1] = 5;
		corner_index[2] = 6;
		corner_index[3] = 7;
		edge_index[0] = 4;
		edge_index[1] = 5;
		edge_index[2] = 6;
		edge_index[3] = 7;
	}
}

void FMFace::printFace( ostream &stream )
{
	stream << "\tFace: " << index << " Corner_Indecies: { " << corner_index[0] << ", " << corner_index[1] << ", " << corner_index[2] << ", " << corner_index[3]  << " }" << endl;
}

/**************************************************************************************************
 *
 *	FlowMapCell
 *		The cell in a given flow map contains 6 faces and a bounding box
 *
 *
**************************************************************************************************/
class FlowMapCell
{
  public:
	double bbox[6];
	FMFace face[6];
	FMCornerPoint corner_points[8];

	void INIT();
	void setBB( double *bb );
	void setCorners();

	int getCellFace( double x, double y, double z, int& fi, int& ei, int& ci );
	int computeAcceptableFaces();
	void computeOUTCellFaces();

	int advectOnFlow( double &x, double &y, double &z, double &t, int &f, int &cc);

	void printCell( ostream &stream );

};

void FlowMapCell::INIT()
{
	for( int i = 0; i < 6; i++ )
	{
		face[i].INIT(i);
	}
}

void FlowMapCell::setBB( double *bb  )
{
	bbox[0] = bb[0];
	bbox[1] = bb[1];
	bbox[2] = bb[2];
	bbox[3] = bb[3];
	bbox[4] = bb[4];
	bbox[5] = bb[5];
}

void FlowMapCell::setCorners()
{
	double coordsX[8] = { bbox[0], bbox[1], bbox[0], bbox[1], bbox[0], bbox[1], bbox[0], bbox[1] };
	double coordsY[8] = { bbox[2], bbox[2], bbox[3], bbox[3], bbox[2], bbox[2], bbox[3], bbox[3] };
	double coordsZ[8] = { bbox[4], bbox[4], bbox[4], bbox[4], bbox[5], bbox[5], bbox[5], bbox[5] };

	for( int  cp = 0; cp < 8; cp++ )
	{
		corner_points[cp].SetINPoint( coordsX[cp], coordsY[cp], coordsZ[cp] );
		corner_points[cp].SetINCornerID( cp );
	}
}

int FlowMapCell::computeAcceptableFaces()
{
	int num_accepted = 0;

	for( int f = 0; f < 6; f++ )
	{

		int faceCount[6] = {0,0,0,0,0,0};

		for( int i = 0; i < 4; i++ )
		{
			int newFaceout = corner_points[ face[f].corner_index[i]].out.face_index;
			int newCorner  = corner_points[ face[f].corner_index[i]].out.corner_index;
			int newEdgeID  = corner_points[ face[f].corner_index[i]].out.edge_index;

			if( newFaceout == -1 )
			{
				if( newEdgeID == -1 )
				{
					if( newCorner == -1 )
					{
						break;
					}
					else
					{
						for( int p = 0; p < 3; p++ )
						{
							faceCount[ PointsToFaces[newCorner][p] ]++;
						}
					}
				}
				else
				{
					for( int e = 0; e < 2; e++)
					{
						faceCount[ EdgesToFaces[newEdgeID][e] ]++;
					}
				}
			}
			else
			{
				faceCount[ newFaceout ]++;
			}

		}
		
//		std::cerr << "FaceCount: " << f << ": [ "<< faceCount[0] << ", " <<  faceCount[1] << ", " <<  faceCount[2] << ", " <<  faceCount[3] << ", " <<  faceCount[4] << ", " <<  faceCount[5] << " ]" << std::endl; 

		for( int i = 0; i < 6; i++)
		{
			if( faceCount[i] == 4 )
			{
				face[i].accepted_flow = true;
				num_accepted++;
				break;
			}
		}
	}

	return num_accepted;
}

int FlowMapCell::getCellFace( double x, double y, double z, int& fi, int& ei, int& ci )
{
	
	FMParticle test_particle;
	test_particle.x = x;
	test_particle.y = y;
	test_particle.z = z;

	bool hitCorner = false;
	for( int cp2 = 0; cp2 < 8; cp2++ )
	{

		double xx = corner_points[cp2].in.x;
		double yy = corner_points[cp2].in.y;
		double zz = corner_points[cp2].in.z;


		if( test_particle == corner_points[cp2].in )
		{
			hitCorner = true;
			test_particle.corner_index = corner_points[cp2].in.corner_index;
			break;
		}
	}	

	if( !hitCorner )
	{

		if(      y == bbox[2] && z == bbox[4] )
		{
			test_particle.edge_index = 0;
		} 
		else if( x == bbox[1] && z == bbox[4] )
		{
			test_particle.edge_index = 1;
		} 
		else if( y == bbox[3] && z == bbox[4] )
		{
			test_particle.edge_index = 2;
		} 
		else if( x == bbox[0] && z == bbox[4] )
		{
			test_particle.edge_index = 3;
		} 
		else if( y == bbox[2] && z == bbox[5] )
		{
			test_particle.edge_index = 4;
		} 
		else if( x == bbox[1] && z == bbox[5] )
		{
			test_particle.edge_index = 5;
		} 
		else if( y == bbox[3] && z == bbox[5] )
		{
			test_particle.edge_index = 6;
		} 
		else if( x == bbox[0] && z == bbox[5] )
		{
			test_particle.edge_index = 7;
		} 
		else if( x == bbox[0] && y == bbox[2] )
		{
			test_particle.edge_index = 8;
		} 
		else if( x == bbox[1] && y == bbox[2] )
		{
			test_particle.edge_index = 9;
		} 
		else if( x == bbox[0] && y == bbox[3] )
		{
			test_particle.edge_index = 10;
		} 
		else if( x == bbox[1] && y == bbox[3] )
		{
			test_particle.edge_index = 11;
		}
		else
		{
			if(      x == bbox[0] )
			{
				test_particle.face_index = 3;
			}
			else if( x == bbox[1] )
			{
				test_particle.face_index = 1;
			}
			else if( y == bbox[2] )
			{
				test_particle.face_index = 0;
			}
			else if( y == bbox[3] )
			{
				test_particle.face_index = 2;
			}
			else if( z == bbox[4] )
			{
				test_particle.face_index = 4;
			}
			else if( z == bbox[5] )
			{
				test_particle.face_index = 5;
			}
			else
			{
				cerr << "Why not: " << x << " " << y << " " << z << std::endl;
			}
		}
	}
	if( test_particle.face_index == -1 )
	{
		if( test_particle.corner_index != -1 )
		{
			ci = test_particle.corner_index;
			return PointsToFaces[ test_particle.corner_index][0];
		}
		else if( test_particle.edge_index != -1 )
		{
			ei = test_particle.edge_index;	
			return EdgesToFaces[ test_particle.edge_index][0];	
		}	
	}
	else
	{
		fi = test_particle.face_index;
		return test_particle.face_index;
	}

	return -1;
}

void FlowMapCell::computeOUTCellFaces()
{
	for( int cp = 0; cp < 8; cp++ )
	{
		double x = corner_points[cp].out.x;
		double y = corner_points[cp].out.y;
		double z = corner_points[cp].out.z;

		int fi = -1;
		int ei = -1;
		int ci = -1;

		int f = getCellFace( x, y, z, fi, ei, ci  );	
		corner_points[cp].out.face_index = fi;
		corner_points[cp].out.edge_index = ei;
		corner_points[cp].out.corner_index = ci;
	}
}
#if 0
		bool hitCorner = false;
		for( int cp2 = 0; cp2 < 8; cp2++ )
		{
			if( corner_points[cp].out == corner_points[cp2].in )
			{
				hitCorner = true;
				corner_points[cp].out.corner_index = corner_points[cp2].in.corner_index;
				break;
			}
		}	

		if( !hitCorner )
		{
			double x = corner_points[cp].out.x;
			double y = corner_points[cp].out.y;
			double z = corner_points[cp].out.z;

			if(      y == bbox[2] && z == bbox[4] )
			{
				corner_points[cp].out.edge_index = 0;
			} 
			else if( x == bbox[1] && z == bbox[4] )
			{
				corner_points[cp].out.edge_index = 1;
			} 
			else if( y == bbox[3] && z == bbox[4] )
			{
				corner_points[cp].out.edge_index = 2;
			} 
			else if( x == bbox[0] && z == bbox[4] )
			{
				corner_points[cp].out.edge_index = 3;
			} 
			else if( y == bbox[2] && z == bbox[5] )
			{
				corner_points[cp].out.edge_index = 4;
			} 
			else if( x == bbox[1] && z == bbox[5] )
			{
				corner_points[cp].out.edge_index = 5;
			} 
			else if( y == bbox[3] && z == bbox[5] )
			{
				corner_points[cp].out.edge_index = 6;
			} 
			else if( x == bbox[0] && z == bbox[5] )
			{
				corner_points[cp].out.edge_index = 7;
			} 
			else if( x == bbox[0] && y == bbox[2] )
			{
				corner_points[cp].out.edge_index = 8;
			} 
			else if( x == bbox[1] && y == bbox[2] )
			{
				corner_points[cp].out.edge_index = 9;
			} 
			else if( x == bbox[0] && y == bbox[3] )
			{
				corner_points[cp].out.edge_index = 10;
			} 
			else if( x == bbox[1] && y == bbox[3] )
			{
				corner_points[cp].out.edge_index = 11;
			}
			else
			{
				if(      x == bbox[0] )
				{
					corner_points[cp].out.face_index = 0;
				}
				else if( x == bbox[1] )
				{
					corner_points[cp].out.face_index = 1;
				}
				else if( y == bbox[1] )
				{
					corner_points[cp].out.face_index = 2;
				}
				else if( y == bbox[1] )
				{
					corner_points[cp].out.face_index = 3;
				}
				else if( z == bbox[1] )
				{
					corner_points[cp].out.face_index = 4;
				}
				else if( z == bbox[1] )
				{
					corner_points[cp].out.face_index = 5;
				}
				else
				{
					cerr << "Why not: " << x << " " << y << " " << z << std::endl;
					
					//cerr << "No Face or edge or corner" << endl;
				}
			}


		}

	}
}

#endif

int FlowMapCell::advectOnFlow( double &x, double &y, double &z, double &t, int &f, int &cc )
{
	if( !face[f].accepted_flow )
	{
		return 2; //advection not supported
	}

/* Find Fractional Weights for interpolation  */
	double fracD1;
	double fracD2;

	FMParticle p1 = corner_points[ face[ f ].corner_index[0] ].in;
	FMParticle p2 = corner_points[ face[ f ].corner_index[1] ].in;
	FMParticle p3 = corner_points[ face[ f ].corner_index[2] ].in;
	FMParticle p4 = corner_points[ face[ f ].corner_index[3] ].in;

	if( p2.x != p1.x  )
	{
		fracD1 = ( x - p1.x) / ( p2.x - p1.x);
	}
	else if( p2.y != p1.y )
	{
		fracD1 = ( y - p1.y) / ( p2.y - p1.y);
	}
	else if( p2.z != p1.z )
	{
		fracD1 = ( z - p1.z) / ( p2.z - p1.z);
	}
	else
	{
		return 2; //Error state should not occur but we can always euler
	}

	if( p3.x != p1.x  )
	{
		fracD2 = ( x - p1.x) / ( p3.x - p1.x);
	}
	else if( p3.y != p1.y )
	{
		fracD2 = ( y - p1.y) / ( p3.y - p1.y);
	}
	else if( p3.z != p1.z )
	{
		fracD2 = ( z - p1.z) / ( p3.z - p1.z);
	}
	else
	{
		return 2; //Error state should not occur but we can always euler
	}

cerr << endl;
cerr << "Particle_Advect:" << endl;
cerr << "P: " << x << " " << y << " " << z << endl;
cerr << endl;
cerr << "Particle_INPUT Points:" << endl;
cerr << "P1: " << p1.x << " " << p1.y << " " << p1.z << endl;
cerr << "P2: " << p2.x << " " << p2.y << " " << p2.z << endl;
cerr << "P3: " << p3.x << " " << p3.y << " " << p3.z << endl;
cerr << "P4: " << p4.x << " " << p4.y << " " << p4.z << endl;
cerr << endl;

	p1 = corner_points[ face[ f ].corner_index[0] ].out;
	p2 = corner_points[ face[ f ].corner_index[1] ].out;
	p3 = corner_points[ face[ f ].corner_index[2] ].out;
	p4 = corner_points[ face[ f ].corner_index[3] ].out;

cerr << "Particle_OUTPUT Points:" << endl;
cerr << "P1: " << p1.x << " " << p1.y << " " << p1.z << endl;
cerr << "P2: " << p2.x << " " << p2.y << " " << p2.z << endl;
cerr << "P3: " << p3.x << " " << p3.y << " " << p3.z << endl;
cerr << "P4: " << p4.x << " " << p4.y << " " << p4.z << endl;
cerr <<endl;
cerr << "FracD1: " << fracD1 << endl;
cerr << "FracD2: " << fracD2 << endl;

	double X01 = p1.x + fracD1 * ( p2.x - p1.x ); 
	double X23 = p1.y + fracD1 * ( p4.x - p3.x );

	double Y01 = p1.y + fracD1 * ( p2.y - p1.y ); 
	double Y23 = p3.y + fracD1 * ( p4.y - p3.y );

	double Z01 = p1.z + fracD1 * ( p2.z - p1.z ); 
	double Z23 = p3.z + fracD1 * ( p4.z - p3.z );

	double T01 = p1.time + fracD1 * ( p2.time - p1.time ); 
	double T23 = p3.time + fracD1 * ( p4.time - p3.time );

	double tt = t;
	double xx = x;
	double yy = y;
	double zz = z;
	double ff = f;
	
	x = X01 + fracD2 * ( X23-X01 );
	y = Y01 + fracD2 * ( Y23-Y01 );
	z = Z01 + fracD2 * ( Z23-Z01 );
	t = T01 + fracD2 * ( T23-T01 );
	f = FaceToFace[f];

	if( tt == t )
	{
		std::cerr << "dt = 0 " << std::endl;
		return 0;
	}
	//cell change in which direction
	//define these directions with chart. This is example though I think not right
	//0 -> 2 = -x ( 1 )
	//2 -> 0 = +x ( 2 )
	//1 -> 3 = 
	if( ( f == 0 && ff == 2 ) ||() )
	{
	}

	return 1;

}

void FlowMapCell::printCell( ostream &stream )
{
	stream << "BBOX: [ " << bbox[0] << ", " << bbox[1] << ", " << bbox[2] << ", " << bbox[3] << ", " << bbox[4] << ", " << bbox[5] << " ]" << endl; 

	stream << "Faces: \n    [" << endl;

	for( int f = 0; f < 6; f++ )
	{
		face[f].printFace( stream );
	}
	for( int c = 0; c < 8; c++ )
	{
		stream << "\tCorner Point: " << c << endl;
		corner_points[c].printCorner( stream );
	}

	stream << "    ]" << endl;

}

#endif
