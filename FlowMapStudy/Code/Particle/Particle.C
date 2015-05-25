#include <Particle.h>

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


double Particle::getStepSize()
{
	return stepSize;
}
void Particle::setStepSize( double ss)
{ 
	stepSize = ss;
}


int getEdgeIDs( int* bb  )
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

int ComputeFID( int* bb )
{

const int BoundsToFace[6] = { 3, 1, 0, 2, 4, 5 };

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
		int edgeID = getEdgeIDs( bb );
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

void Particle::setFID( double* bbox )
{
	int bb[6] = {0,0,0,0,0,0};

	bb[0] = ( ( almostEqual(x, bbox[0] )) ? 1 : 0 );
	bb[1] = ( ( almostEqual(x, bbox[1] )) ? 1 : 0 );
	bb[2] = ( ( almostEqual(y, bbox[2] )) ? 1 : 0 );
	bb[3] = ( ( almostEqual(y, bbox[3] )) ? 1 : 0 );
	bb[4] = ( ( almostEqual(z, bbox[4] )) ? 1 : 0 );
	bb[5] = ( ( almostEqual(z, bbox[5] )) ? 1 : 0 );

	FID = ComputeFID( bb );

	if( FID > 27 )
	{
		cerr << "SetFID: Error: Conditions not met to set a face for flow"  << bb[0] << " " << bb[1] << " " << bb[2] << " " << bb[3] << " " << bb[4] << " " << bb[5] << endl;
		cerr << FID << " " << x << " " << y << " " << z << endl;	
	}
}
