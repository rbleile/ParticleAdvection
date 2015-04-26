#include "Flow.h"
#include "Point.h"

/***************************************************************************************
*
*  All Possible Cell mappings between points <-> faces and edges <-> faces within a cell
*
***************************************************************************************/

//Definition: Which faces belong to each point ( i.e. given a point you know what the connecting faces are in your cell )
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

//Definition: Which points belong to eahc face ( i.e. given a face you know what the points are in your cell )
const int FacesToPoints[6][4] = 
	{
		{0,1,4,5}, //0
		{1,3,5,7}, //1
		{2,3,6,7}, //2
		{0,2,4,6}, //3
		{0,1,2,3}, //4
		{4,5,6,7}  //5
	};

//Definition: Which faces belong to which edges ( i.e. given an edge you know what the connceting faces are in your cell )
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
		{2, 3}
	};

//Definition: Which edges belong to which faces ( i.e. given an face you know what the edges are in your cell )
const int FacesToEdges[6][4] = 
	{
		{0, 4, 8,9 },
		{1, 5, 9,10},
		{2, 6,10,11},
		{3, 7, 8,11},
		{0, 1, 2, 3},
		{4, 5, 6, 7}
	};

//Definition: Which face connects two cells ( i.e. given a face in my cell what is the face in a connecting cell )
const int FaceToFace[6]= { 2, 3, 0, 1, 5, 4 };

//Definition: which face lies on which boundry ( xmin, xmax, ymin, ymax, zmin, zmax )
const int BoundsToFace[6] = { 3, 1, 0, 2, 4, 5 };

inline bool operator==( const Point &p1, const Point &p2 )
{
	if(
	      fabs((p1.x-p2.x)) <= epsilon &&
		  fabs((p1.y-p2.y)) <= epsilon &&
		  fabs((p1.z-p2.z)) <= epsilon
      )
	{
		return true;
	}
	else
	{
		return false;
	}
}

inline bool operator!=( const Point &p1, const Point &p2 )
{
	return !(p1==p2);
}

//Function for checking floating point equvalence to epsilon error
inline bool almostEqual( float p1, float p2 )
{
	return ( ( fabs((p1-p2)) <= epsilon ) ? true : false);
} 

inline bool almostEqual( double p1, double p2 )
{
	return ( ( fabs((p1-p2)) <= epsilon ) ? true : false);
} 

int getEdgeID( int* bb  )
{

	if     ( bb[2] && bb[4] ) return 0;	
	else if( bb[1] && bb[4] ) return 1;	
	else if( bb[3] && bb[4] ) return 2;	
	else if( bb[0] && bb[4] ) return 3;	
	else if( bb[2] && bb[5] ) return 4;	
	else if( bb[1] && bb[5] ) return 5;	
	else if( bb[3] && bb[5] ) return 6;	
	else if( bb[0] && bb[5] ) return 7;	
	else if( bb[0] && bb[2] ) return 8;	
	else if( bb[1] && bb[2] ) return 9;	
	else if( bb[1] && bb[3] ) return 10;
	else if( bb[0] && bb[3] ) return 11;
	else cerr << "getEdgeID: Error: No valid edge detected." << endl;

	return -15;

}

int computeFID( int* bb )
{
	int count = 0;
	for( int i = 0; i < 6; i++ )
	{
		count += bb[i];
	}

	if( count == 0 )
	{
		return -1;	
	}
	else if( count == 1 )
	{
		for( int i = 0; i < 6; i++ )
		{
			if( bb[i] == 1  )
			{
				return BoundsToFace[i];
			}
		}
	}
	else if( count == 2 )
	{
		int edgeID = getEdgeID( bb );
		return 14+edgeID;
	}
	else if( count == 3 )
	{
		int pointID = 0;
		pointID += ( ( bb[0] ) ? 1 : 0 );
		pointID += ( ( bb[2] ) ? 2 : 0 );
		pointID += ( ( bb[4] ) ? 4 : 0 );

		return 6+pointID;
	}
	else
	{
		cerr << "computeFID: Error: cannot have this case. Touching too many boundaries" << endl;
	}

	return -1;
}

void Flow::printFlow( ostream &stream )
{
	stream << "CellID: " << cellID << "\t";
	stream << "FID: " << FID << "\t";
	stream << "In Point: < X: " << in.x << ", Y: " << in.y << ", Z: " << in.z << ", T: " << in.t << " > " << "\t";
	stream << "Out Point: < X: " << out.x << ", Y: " << out.y << ", Z: " << out.z << ", T: " << out.t << " > " << endl;
}
void Flow::setPositions( double x, double y, double z )
{
	in.setPoint( x, y, z );
	out.setPoint( x, y, z );
}
void Flow::setFID( double* bbox )
{
	int bb[6] = {0,0,0,0,0,0};

	bb[0] = ( ( almostEqual(out.x, bbox[0] )) ? 1 : 0 );
	bb[1] = ( ( almostEqual(out.x, bbox[1] )) ? 1 : 0 );
	bb[2] = ( ( almostEqual(out.y, bbox[2] )) ? 1 : 0 );
	bb[3] = ( ( almostEqual(out.y, bbox[3] )) ? 1 : 0 );
	bb[4] = ( ( almostEqual(out.z, bbox[4] )) ? 1 : 0 );
	bb[5] = ( ( almostEqual(out.z, bbox[5] )) ? 1 : 0 );

	FID = ( out == in ) ? -1 : computeFID( bb );

	if( FID > 27 )
	{
		cerr << "SetFID: Error: Conditions not met to set a face for flow"  << bb[0] << " " << bb[1] << " " << bb[2] << " " << bb[3] << " " << bb[4] << " " << bb[5] << endl;
		cerr << FID << " " << out.x << " " << out.y << " " << out.z << endl;	
	}

/*
	if( FID == -1 )
	{
		cerr << "SetFID: Error: Conditions not met to set a face for flow"  << bb[0] << " " << bb[1] << " " << bb[2] << " " << bb[3] << " " << bb[4] << " " << bb[5] << endl;
		cerr << FID << " " << out.x << " " << out.y << " " << out.z << endl;
	}
*/
}


