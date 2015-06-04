#include <DomainMesh.h>

void DomainMesh::fillInFlows()
{

	for( int d = 0; d < 3; d++ )
	{
		Flow* flows = flowField[d];

		long int Rdims[3];
		if( d == 0 )     { Rdims[0] = nx;		Rdims[1] = mesh->ny;	Rdims[2] = mesh->nz; }
		else if( d == 1 ){ Rdims[0] = mesh->nx; Rdims[1] = ny;			Rdims[2] = mesh->nz; }
		else if( d == 2 ){ Rdims[0] = mesh->nx; Rdims[1] = mesh->ny;	Rdims[2] = nz; }

		double Rpos[3] = { x0, y0, z0 };

		double Rdel[3];
		if( d == 0 )     { Rdel[0] = dx;		Rdel[1] = mesh->dy; Rdel[2] = mesh->dz; }
		else if( d == 1 ){ Rdel[0] = mesh->dx;	Rdel[1] = dy;		Rdel[2] = mesh->dz; }
		else if( d == 2 ){ Rdel[0] = mesh->dx;	Rdel[1] = mesh->dy; Rdel[2] = dz; }

		for( long int z_id = 0; z_id < Rdims[2]; z_id++ )
		{	
			for( long int y_id = 0; y_id < Rdims[1]; y_id++ )
			{	
				for( long int x_id = 0; x_id < Rdims[0]; x_id++ )
				{	
	
					long int id = x_id + Rdims[0] * ( y_id + Rdims[1] * z_id );

					double x = Rpos[0] + Rdel[0]*x_id;
					double y = Rpos[1] + Rdel[1]*y_id;
					double z = Rpos[2] + Rdel[2]*z_id;

					flows[id].setPositions( x, y, z );
				}
			}
		}
	}
}

void DomainMesh::setFlowsCellIDs()
{
	double xmax = x0 + dx*(nx-1);
	double ymax = y0 + dy*(ny-1);
	double zmax = z0 + dz*(nz-1);

	for( int d = 0; d < 3; d++ )
	{
		Flow* flows = flowField[d];

		long int Rdims[3];
		if( d == 0 )     { Rdims[0] = nx;		Rdims[1] = mesh->ny; Rdims[2] = mesh->nz; }
		else if( d == 1 ){ Rdims[0] = mesh->nx; Rdims[1] = ny;		 Rdims[2] = mesh->nz; }
		else if( d == 2 ){ Rdims[0] = mesh->nx; Rdims[1] = mesh->ny; Rdims[2] = nz; }

		for( long int z_id = 0; z_id < Rdims[2]; z_id++ )
		{	
			for( long int y_id = 0; y_id < Rdims[1]; y_id++ )
			{	
				for( long int x_id = 0; x_id < Rdims[0]; x_id++ )
				{	
	
					long int id = x_id + Rdims[0] * ( y_id + Rdims[1] * z_id );

					Point tp = flows[id].in;

					double dt = 0.001;
					double vel[3];
					mesh->getVelocity( tp, vel );

					Point np( tp.x + vel[0]*dt, tp.y + vel[1]*dt, tp.z + vel[2]*dt );	

					if( np.x > xmax || np.x < x0 
					 || np.y > ymax || np.y < y0
					 || np.z > zmax || np.z < z0 )
					{
						flows[id].cellID = -1;
					}
					else
					{
						flows[id].cellID = getCellID( np.x, np.y, np.z );

						if( np.x == xmax  )
						{
							flows[id].cellID -= 1;
						}
						if( np.y == ymax  )
						{
							flows[id].cellID -= nx;
						}
						if( np.z == zmax  )
						{
							flows[id].cellID -= nx*ny;
						}		
					}
				}
			}
		}
	}
}

#include <cstdio>
#include <iomanip>

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

int DomainMesh::EulerCellAdvection( long int cellID, double endTime, double* bb, Particle &particle )
{
	int stat = 1;

	double cell_bb[6];
	getCellBounds( cellID, cell_bb );

	while( stat == 1 )
	{
		stat = mesh->RK4( cell_bb, bb, endTime, particle );
	}

	return stat;
} 

int DomainMesh::ReverseEulerCellAdvection( long int cellID, double endTime, double* bb, Particle &particle )
{
	int stat = 1;

	double cell_bb[6];
	getCellBounds( cellID, cell_bb );

	while( stat == 1 )
	{
		//stat = mesh->REV_Euler( cell_bb, bb, endTime, particle );
		stat = mesh->REV_RK4( cell_bb, bb, endTime, particle );
	}

	return 0;
} 

long int DomainMesh::computeAllAcceptableFlows()
{
	long int acceptableFlows = 0;
	long int no_count_face = 0;
	long int total = 0;

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


	for( int d = 0; d < 3; d++ )
	{

		int d_id1 = (d + 1)%3;
		int d_id2 = (d + 2)%3;

		Flow* flows = flowField[d];

		long int Rdims[3];
		if( d == 0 )     { Rdims[0] = nx;		Rdims[1] = mesh->ny;	Rdims[2] = mesh->nz; }
		else if( d == 1 ){ Rdims[0] = mesh->nx;	Rdims[1] = ny;			Rdims[2] = mesh->nz; }
		else if( d == 2 ){ Rdims[0] = mesh->nx;	Rdims[1] = mesh->ny;	Rdims[2] = nz; }

		long int Rcdim[3];
		if( d == 0 )     { Rcdim[0] = nx;			Rcdim[1] = mesh->ny-1;	Rcdim[2] = mesh->nz-1; }
		else if( d == 1 ){ Rcdim[0] = mesh->nx-1;	Rcdim[1] = ny;			Rcdim[2] = mesh->nz-1; }
		else if( d == 2 ){ Rcdim[0] = mesh->nx-1;	Rcdim[1] = mesh->ny-1;	Rcdim[2] = nz; }

		long int xstep = 1;
		long int ystep = Rdims[0];
		long int zstep = Rdims[0]*Rdims[1];

		long int Rneighbors[3];
		if( d == 0 )     { Rneighbors[0] = ystep;	Rneighbors[1] = zstep;	Rneighbors[2] = ystep+zstep; }
		else if( d == 1 ){ Rneighbors[0] = zstep;	Rneighbors[1] = xstep;	Rneighbors[2] = zstep+xstep; }
		else if( d == 2 ){ Rneighbors[0] = xstep;	Rneighbors[1] = ystep;	Rneighbors[2] = xstep+ystep; }

		for( long int z_id = 0; z_id < Rcdim[2]; z_id++ )
		{	
			for( long int y_id = 0; y_id < Rcdim[1]; y_id++ )
			{	
				for( long int x_id = 0; x_id < Rcdim[0]; x_id++ )
				{

					long int cell_id = x_id + (Rdims[0]-1) * ( y_id + (Rdims[1]-1) * z_id );
					long int base_flow = x_id + Rdims[0] * ( y_id + Rdims[1] * z_id );

					long int flowIDs[4] = { base_flow, base_flow + Rneighbors[0], base_flow + Rneighbors[1], base_flow + Rneighbors[2] };
					bool isGood = true;
					int l_FID[3] = {-1, -1, -1};
					long int matchCell = flows[ flowIDs[0] ].cellID;

					total++;

					int canCompute = 4;
					for( int p = 0; p < 4; p++ )
					{
						double vel[3];
						mesh->getVelocity( flows[ flowIDs[p] ].in, vel );
						if( vel[0] == 0 && vel[1] == 0 && vel[2] == 0 )
						{
							canCompute--;
						}
					}
					if( canCompute==0 )
					{
						no_count_face++;
						AcceptableFlow[d][ cell_id ] = false;
						continue;	
					}

					for( int p = 0; p < 4; p++ )
					{
						int fid = flows[ flowIDs[p] ].FID;
						int cid = flows[ flowIDs[p] ].cellID;

						if( cid != matchCell ){ isGood = false; break; }
						if( fid < 0 ){ isGood = false; break; }
						if( p == 0 )
						{
							if( fid < 6 )
							{
								l_FID[0] = fid;
							}
							else if( fid > 14 )
							{
								int edgeID = fid-14;
								l_FID[0] = EdgesToFaces[edgeID][0];
								l_FID[1] = EdgesToFaces[edgeID][1];
							}
							else
							{
								int pointID = fid - 6;
								l_FID[0] = PointsToFaces[pointID][0];
								l_FID[1] = PointsToFaces[pointID][1];
								l_FID[2] = PointsToFaces[pointID][2];
							}

						}
						else
						{
							if( fid < 6 )
							{
								for( int i = 0; i < 3; i++ )
								{
									if( l_FID[i] != fid )
									{
										l_FID[i] = -1;
									}
								}
							}
							else if( fid > 14 )
							{
								int edgeID = fid-14;
								for( int i = 0; i < 3; i++ )
								{
									if( l_FID[i] == -1) continue;
									bool match = false;
									for( int j = 0; j < 2; j++ )
									{
										if( l_FID[i] == EdgesToFaces[edgeID][j] )
										{
											match = true;
										}
									}
									if( !match )
									{
										l_FID[i] = -1;
									}
								}
							}
							else
							{
								int pointID = fid-6;
								for( int i = 0; i < 3; i++ )
								{
									if( l_FID[i] == -1) continue;
									bool match = false;
									for( int j = 0; j < 3; j++ )
									{
										if( l_FID[i] == PointsToFaces[pointID][j] )
										{
											match = true;
										}
									}
									if( !match )
									{
										l_FID[i] = -1;
									}
								}
							}
						}
					}
					
					if( isGood )
					{
						bool isComplete = false;
						for( int i = 0; i < 3; i++ )
						{
							if( l_FID[i] != -1 )
							{
								isComplete = true;
							}
						}
						isGood = isComplete;
					}

					AcceptableFlow[d][ cell_id ] = isGood;
					acceptableFlows += ( isGood ) ? 1 : 0;
				}
			}
		}
	}

	cerr << "Throw Away: " << no_count_face << endl;
	cerr << "Old Total Num Faces: " << total << endl; 
	cerr << "New Total Num Faces: " << total-no_count_face << endl; 
	cerr << "Num Computable Faces: " << acceptableFlows << endl;
	cerr << "Percent Accepted New: " << ((double)acceptableFlows) / ((double)(total-no_count_face)) * 100 << "%" << endl; 

	return acceptableFlows;
}

