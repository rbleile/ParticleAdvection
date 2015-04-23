
#include "FlowMap.h"

#ifndef MESH_H
#define MESH_H

//*************************************************************************************************
//    Mesh Class for Rectilinear or Regular Meshes
//*************************************************************************************************
class Mesh
{
  public:

    // Is the Mesh Regular or Rectilinear
    bool isRegular;

    // Number of Points in Each Dimension
    int nx;
    int ny;
    int nz;

    // Step Size in Each Dimension
    double dx;
    double dy;
    double dz;

    // Origin Point of Mesh
    double x0;
    double y0;
    double z0;

    // Mesh Vector Field
    double* V;

    // Mesh Fiels Points ( Rectilinear Representation )
    double* X;
    double* Y;
    double* Z;

	// Flow Map Cells ( used for flow maps )
	int num_cells;
	FlowMapCell* fmc;

    //Mesh Empty Mesh Constructor
    Mesh() : isRegular( false ), X(NULL), Y(NULL), Z(NULL), V(NULL), nx(0), ny(0), nz(0), dx(-1), dy(-1), dz(-1), x0(0), y0(0), z0(0){};

    //Mesh Regular Mesh Constructor
    Mesh( bool isReg, int num_x, int num_y, int num_z, double* v_field, double del_x, double del_y, double del_z, double x_p, double y_p, double z_p ) 
		: isRegular( isReg ), X(NULL), Y(NULL), Z(NULL), V(v_field), nx(num_x), ny(num_y), nz(num_z), dx(del_x), dy(del_y), dz(del_z), x0(x_p), y0(y_p), z0(z_p){};

    //Mesh Rectilinear Mesh Constructor
    Mesh( bool isReg, int num_x, int num_y, int num_z, double* x_points, double* y_points, double* z_points, double* v_field, double x_p, double y_p, double z_p ) 
		: isRegular( isReg ), X(x_points), Y(y_points), Z(z_points), V(v_field), nx(num_x), ny(num_y), nz(num_z), dx(-1), dy(-1), dz(-1), x0(x_p), y0(y_p), z0(z_p){};

    //Mesh Full Mesh Constructor
    Mesh( bool isReg, int num_x, int num_y, int num_z, double* x_points, double* y_points, double* z_points, double* v_field, double del_x, double del_y, double del_z, double x_p, double y_p, double z_p ) 
		: isRegular( isReg ), X(x_points), Y(y_points), Z(z_points), V(v_field), nx(num_x), ny(num_y), nz(num_z), dx(del_x), dy(del_y), dz(del_z), x0(x_p), y0(y_p), z0(z_p){};

     // Get the Cell ID for a given Point in a given coordinate
     int getRegularCellID_X( double xp ){
         return (int)((xp-x0)/dx);
     }

     int getRegularCellID_Y( double yp ){
         return (int)((yp-y0)/dy);
     }

     int getRegularCellID_Z( double zp ){
         return (int)((zp-z0)/dz);
     }

     void getRegularLogicalCellID( double xp, double yp, double zp, int *cellID )
     {
         cellID[0] = getRegularCellID_X( xp );
         cellID[1] = getRegularCellID_Y( yp );
         cellID[2] = getRegularCellID_Z( zp );
     }

     // Simple Rectilinear Cell ID LookUp and Calculation 
     void getRectilinearLogicalCellID( double xp, double yp, double zp, int *cellID )
     {
		cellID[0] = 0;
		cellID[1] = 0;
		cellID[2] = 0;

         for( int i = 1; i < nx; i++ )
         {
             if( xp < X[i] )
             {
                 cellID[0] = i-1;
             }
         }

         for( int j = 1; j < ny; j++ )
         {
             if( yp < Y[j] )
             {
                 cellID[1] = j-1;
             }
         }

         for( int k = 1; k < nz; k++ )
         {
             if( zp < X[k] )
             {
                 cellID[2] = k-1;
             }
         }

     }
     void getLogicalCellID( double xp, double yp, double zp, int *cellID )
     {
 
       if(isRegular)
         {
             getRegularLogicalCellID( xp, yp, zp, cellID );
         }
		else{
             getRectilinearLogicalCellID( xp, yp, zp, cellID );
		}
     }

     // Simple Regular Cell ID Calculation
	int getRegularCellID( double xp, double yp, double zp )
    {
		int cellID[3];

		getRegularLogicalCellID( xp, yp, zp, cellID );

        return cellID[2]*(nx-1)*(ny-1)+cellID[2]*(nx-1)+cellID[0];
    }

     // Simple Rectilinear Cell ID LookUp and Calculation 
     int getRectilinearCellID( double xp, double yp, double zp )
     {
		int cellID[3];

		getRectilinearLogicalCellID( xp, yp, zp, cellID );

        return cellID[2]*(nx-1)*(ny-1)+cellID[2]*(nx-1)+cellID[0];
     }

     //Generic Function to compute the cell id of a point on this mesh
     int getCellID( double xp, double yp, double zp )
     {
         if( isRegular )
             return getRegularCellID( xp, yp, zp );
         else
           return getRectilinearCellID( xp, yp, zp );
     }

	void getCellXRangeRegular( int xid, double &xmin, double &xmax )
	{
		xmin = x0 + dx*xid;
		xmax = x0 + dx*(xid+1);	
	}
	void getCellYRangeRegular( int yid, double &ymin, double &ymax )
	{
		ymin = y0 + dy*yid;
		ymax = y0 + dy*(yid+1);	
	}
	void getCellZRangeRegular( int zid, double &zmin, double &zmax )
	{
		zmin = z0 + dz*zid;
		zmax = z0 + dz*(zid+1);	
	}

	void getCellXRangeRectilinear( int xid, double &xmin, double &xmax )
	{
		xmin = x0 + dx*xid;
		xmax = x0 + dx*(xid+1);	
	}
	void getCellYRangeRectilinear( int yid, double &ymin, double &ymax )
	{
		ymin = y0 + dy*yid;
		ymax = y0 + dy*(yid+1);	
	}
	void getCellZRangeRectilinear( int zid, double &zmin, double &zmax )
	{
		zmin = z0 + dz*zid;
		zmax = z0 + dz*(zid+1);	
	}
	//Get the Min and Max of a cell given a logical id
	void getCellXRange( int xid, double &xmin, double &xmax )
	{
		if( isRegular )
			getCellXRangeRegular( xid, xmin, xmax );
		else
			getCellXRangeRectilinear( xid, xmin, xmax );
	}
	void getCellYRange( int yid, double &ymin, double &ymax )
	{
		if( isRegular )
			getCellYRangeRegular( yid, ymin, ymax );
		else
			getCellYRangeRectilinear( yid, ymin, ymax );
	}
	void getCellZRange( int zid, double &zmin, double &zmax )
	{
		if( isRegular )
			getCellZRangeRegular( zid, zmin, zmax );
		else
			getCellZRangeRectilinear( zid, zmin, zmax );
	}

};

void PrintFlowMap( Mesh *mesh, ostream &stream )
{
	for( int i = 0; i < mesh->num_cells; i++ )
	{
		stream << "Cell: " << i << endl;
		mesh->fmc[i].printCell( stream );
	}
}

void fillCellParticlesIN( Mesh* mesh  )
{

	int npx = mesh->nx;
	int npy = mesh->ny;
	int npz = mesh->nz;

	int ncx = npx - 1;
	int ncy = npy - 1;
	int ncz = npz - 1;

	for( int k = 0; k < ncz; k++ ){
		for( int j = 0; j < ncy; j++ ){
			for( int i = 0; i < ncx; i++ ){

				int cellID = i + ncx*j + ncx*ncy*k;

				FlowMapCell &Lfmc = mesh->fmc[cellID];

				Lfmc.INIT();

				double xmin, xmax;				
				double ymin, ymax;				
				double zmin, zmax;				

				mesh->getCellXRange( i, xmin, xmax );
				mesh->getCellYRange( j, ymin, ymax );
				mesh->getCellZRange( k, zmin, zmax );

				double BBOX[6] = { xmin, xmax, ymin, ymax, zmin, zmax };

				Lfmc.setBB( BBOX );

				Lfmc.setCorners();

			}
		}
	}
}

#endif
