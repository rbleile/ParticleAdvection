#ifndef DOMAINMESH_H
#define DOMAINMESH_H

#include "Mesh.h"
#include "Flow.h"
#include "Particle.h"

class DomainMesh
{
  public:

    // Number of Points in Each Dimension
    long int nx;
    long int ny;
    long int nz;

    // Step Size in Each Dimension
    double dx;
    double dy;
    double dz;

    // Origin Point of Mesh
    double x0;
    double y0;
    double z0;

	// Flows
	Flow* flowField[3];
	bool* AcceptableFlow[3];
	double* max_diff[3];

	//Fine Grain Mesh
	Mesh* mesh;

    //Mesh Empty Mesh Constructor
    DomainMesh() : nx(0), ny(0), nz(0), dx(0), dy(0), dz(0), x0(0), y0(0), z0(0){};

    //Mesh Regular Mesh Constructor
    DomainMesh( long int num_x, long int num_y, long int num_z, double del_x, double del_y, double del_z, double x_p, double y_p, double z_p, Mesh* fineMesh ) 
		: nx(num_x), ny(num_y), nz(num_z), dx(del_x), dy(del_y), dz(del_z), x0(x_p), y0(y_p), z0(z_p), mesh( fineMesh ){};

     // Get the Cell ID for a given Point in a given coordinate
     long int getRegularCellID_X( double xp ){
         return (int)((xp-x0)/dx);
     }

     long int getRegularCellID_Y( double yp ){
         return (int)((yp-y0)/dy);
     }

     long int getRegularCellID_Z( double zp ){
         return (int)((zp-z0)/dz);
     }

     void getRegularLogicalCellID( double xp, double yp, double zp, long int *cellID )
     {
         cellID[0] = getRegularCellID_X( xp );
         cellID[1] = getRegularCellID_Y( yp );
         cellID[2] = getRegularCellID_Z( zp );
     }

     void getLogicalCellID( double xp, double yp, double zp, long int *cellID )
     {
		getRegularLogicalCellID( xp, yp, zp, cellID );
     }

     // Simple Regular Cell ID Calculation
	long int getRegularCellID( double xp, double yp, double zp )
    {
		long int cellID[3];

		getRegularLogicalCellID( xp, yp, zp, cellID );

		return cellID[0] + (nx-1) * ( cellID[1] + (ny-1) * cellID[2]);
    }

   //Generic Function to compute the cell id of a point on this mesh
   long int getCellID( double xp, double yp, double zp )
   {
	return getRegularCellID( xp, yp, zp );
   }

	void getCellXRange( long int xid, double &xmin, double &xmax )
	{
		xmin = x0 + dx*xid;
		xmax = x0 + dx*(xid+1);	
	}
	void getCellYRange( long int yid, double &ymin, double &ymax )
	{
		ymin = y0 + dy*yid;
		ymax = y0 + dy*(yid+1);	
	}
	void getCellZRange( long int zid, double &zmin, double &zmax )
	{
		zmin = z0 + dz*zid;
		zmax = z0 + dz*(zid+1);	
	}

	//Get the boundingBox for a cell given global id
	void getCellBounds( long int cid, double* bb )
	{
		long int ids[3];
		D1to3C( cid, ids );
		getCellXRange( ids[0], bb[0], bb[1] );
		getCellYRange( ids[1], bb[2], bb[3] );
		getCellZRange( ids[2], bb[4], bb[5] );
	}

	//Convert an index from a one dimensional representaiton to a three dimensional representation
	void D1to3C( long int id, long int *ids )
	{
		ids[0] = id % (nx-1);
		ids[1] = (id/(nx-1)) % (ny-1);
		ids[2] = id / ((nx-1) * (ny-1));
	}

	//Convert an index from a three dimensional representaiton to a one dimensional representation
	long int D3to1C( long int x, long int y, long int z )
	{
		return x + (nx-1)*(y + (ny-1)*z);	
	}

	//Convert an index from a one dimensional representaiton to a three dimensional representation
	void D1to3P( long int id, long int *ids )
	{
		ids[0] = id % nx;
		ids[1] = (id/nx) % ny;
		ids[2] = id / (nx * ny);
	}

	//Convert an index from a three dimensional representaiton to a one dimensional representation
	long int D3to1P( long int x, long int y, long int z )
	{
		return x + nx*(y + ny*z);	
	}

	//Flow functions
	void fillInFlows();
	void setFlowsCellIDs();

	int EulerCellAdvection( long int cellID, double endTime, double* bb, Particle &particle );
	int ReverseEulerCellAdvection( long int cellID, double endTime, double* bb, Particle &particle );

	long int computeAllAcceptableFlows();

};

#endif
