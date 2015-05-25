#include "PrintVTK.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cfloat>
#include <cmath>

using namespace std;

void printVtkFTLE( int nx, int ny, int nz, double xmin, double ymin, double zmin, double dx, double dy, double dz, double *FTLE1, double *FTLE2, double *FTLEDiff )
{

	/*The Uber Mesh*/
	stringstream smesh;
	smesh << "# vtk DataFile Version 3.0" << endl;
	smesh << "vtk output" << endl;
	smesh << "ASCII" << endl;
	smesh << "DATASET RECTILINEAR_GRID" << endl;
	smesh << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
	smesh << "X_COORDINATES " << nx << " float" << endl;
	for( long int i = 0; i < nx; i++ )	smesh << xmin + dx*i << " ";
	smesh << endl;

	smesh << "Y_COORDINATES " << ny << " float" << endl;
	for( long int i = 0; i < ny; i++ )	smesh << ymin + dy*i << " ";
	smesh << endl;

	smesh << "Z_COORDINATES " << nz << " float" << endl;
	for( long int i = 0; i < nz; i++ )	smesh << zmin + dz*i << " ";
	smesh << endl;

	smesh << "POINT_DATA " << nx*ny*nz << endl;
	smesh << "SCALARS FTLEDiff float" << endl;
	smesh << "LOOKUP_TABLE default" << endl;

	double total_diff = 0.0;
	int num_nonz = 0;

	for( int i = 0; i < nx*ny*nz; i++ )
	{
		if( i%9 == 0 && i != 0 )
		{
			smesh << endl;
		}

		total_diff += FTLEDiff[i];
		if( FTLEDiff[i] != 0.0 )
		{
			num_nonz++;
		} 

		smesh << FTLEDiff[i] << " ";

	}

	cerr << "FTLE Non Zero Average Diff: " << total_diff / ( ( num_nonz == 0 ) ? 1 : num_nonz  ) << endl;
	cerr << "FTLE Total Average Diff:    " << total_diff / ( nx*ny*nz) << endl;

	smesh << endl;

	smesh << "FIELD FieldData 3" << endl; 
	smesh << "FTLE1 1 " << nx*ny*nz << " float" << endl;

	for( int i = 0; i < nx*ny*nz; i++ )
	{
		if( i%9 == 0 && i != 0 )
		{
			smesh << endl;
		}
		
		smesh << FTLE1[i] << " ";

	}

	smesh << endl;

	smesh << "FTLE2 1 " << nx*ny*nz << " float" << endl;

	for( int i = 0; i < nx*ny*nz; i++ )
	{
		if( i%9 == 0 && i != 0 )
		{
			smesh << endl;
		}
		
		smesh << FTLE2[i] << " ";

	}

	smesh << endl;

	double total_pdiff = 0.0;
	int num_nonz_p = 0;

	smesh << "FTLE_percent_diff 1 " << nx*ny*nz << " float" << endl;

	for( int i = 0; i < nx*ny*nz; i++ )
	{
		if( i%9 == 0 && i != 0 )
		{
			smesh << endl;
		}

		double avgFTLE = ( FTLE1[i] + FTLE2[i] ) / 2.0;
		double pdiff = (( fabs( FTLE1[i] - FTLE2[i] )) / ((avgFTLE==0) ? 1 : avgFTLE ))*100;

		total_pdiff += pdiff;
		if( pdiff != 0.0 )
		{
			num_nonz_p++;
		} 

		
		smesh << pdiff << " ";

	}

	smesh << endl;

	cerr << "FTLE Non Zero Average Percent Diff: " << total_pdiff / ( ( num_nonz_p == 0 ) ? 1 : num_nonz_p  )  << "%" << endl;
	cerr << "FTLE Total Average Percent Diff:    " << total_pdiff / ( nx*ny*nz) << "%" << endl;

	string fname = "FTLE_OUT.vtk";
	ofstream outM( fname.c_str() );
	outM << smesh.rdbuf();
	outM.close();

}

void VTKFilePrinter::printVtkDs()
{

	/*The Uber Mesh*/
	stringstream smesh;
	smesh << "# vtk DataFile Version 3.0" << endl;
	smesh << "vtk output" << endl;
	smesh << "ASCII" << endl;
	smesh << "DATASET RECTILINEAR_GRID" << endl;
	smesh << "DIMENSIONS " << Udims[0] << " " << Udims[1] << " " << Udims[2] << endl;
	smesh << "X_COORDINATES " << Udims[0] << " float" << endl;
	for( long int i = 0; i < Udims[0]; i++ )	smesh << origin[0] + Udel[0]*i << " ";
	smesh << endl;

	smesh << "Y_COORDINATES " << Udims[1] << " float" << endl;
	for( long int i = 0; i < Udims[1]; i++ )	smesh << origin[1] + Udel[1]*i << " ";
	smesh << endl;

	smesh << "Z_COORDINATES " << Udims[2] << " float" << endl;
	for( long int i = 0; i < Udims[2]; i++ )	smesh << origin[2] + Udel[2]*i << " ";
	smesh << endl;

	stringstream Uname;
	Uname << baseName << "_vtkOut_Mesh.vtk";
	
	string uname;
	Uname >> uname;
	ofstream outM( uname.c_str() );
	outM << smesh.rdbuf();
	outM.close();

	double totAVG = 0;
	int totSum = 0;

	for( int d = 0; d < 3; d++ )
	{
		long int id_1 = (d+1)%3;
		long int id_2 = (d+2)%3;

		long int Rdims[3];
		if     ( d == 0 ){ Rdims[0] = Udims[0]; Rdims[1] = Fdims[1]; Rdims[2] = Fdims[2]; }
		else if( d == 1 ){ Rdims[0] = Fdims[0]; Rdims[1] = Udims[1]; Rdims[2] = Fdims[2]; }
		else if( d == 2 ){ Rdims[0] = Fdims[0]; Rdims[1] = Fdims[1]; Rdims[2] = Udims[2]; }

		long int RMDim[3];
		if     ( d == 0 ){ RMDim[0] = 1;        RMDim[1] = Fdims[1]; RMDim[2] = Fdims[2]; }
		else if( d == 1 ){ RMDim[0] = Fdims[0]; RMDim[1] = 1;        RMDim[2] = Fdims[2]; }
		else if( d == 2 ){ RMDim[0] = Fdims[0]; RMDim[1] = Fdims[1]; RMDim[2] = 1;        }

		long int RcMDim[3];
		if     ( d == 0 ){ RcMDim[0] = 1;			RcMDim[1] = Fdims[1]-1; RcMDim[2] = Fdims[2]-1; }
		else if( d == 1 ){ RcMDim[0] = Fdims[0]-1;	RcMDim[1] = 1;			RcMDim[2] = Fdims[2]-1; }
		else if( d == 2 ){ RcMDim[0] = Fdims[0]-1;	RcMDim[1] = Fdims[1]-1; RcMDim[2] = 1;        }


		double Rdel[3];
		if     ( d == 0 ){ Rdel[0] = Udel[0]; Rdel[1] = Fdel[1]; Rdel[2] = Fdel[2]; }
		else if( d == 1 ){ Rdel[0] = Fdel[0]; Rdel[1] = Udel[1]; Rdel[2] = Fdel[2]; }
		else if( d == 2 ){ Rdel[0] = Fdel[0]; Rdel[1] = Fdel[1]; Rdel[2] = Udel[2]; }


		for( long int uber = 0; uber < Rdims[d]; uber++ )
		{
			stringstream ss;
			ss << "# vtk DataFile Version 3.0" << endl;
			ss << "vtk output" << endl;
			ss << "ASCII" << endl;
			ss << "DATASET RECTILINEAR_GRID" << endl;
			ss << "DIMENSIONS " << RMDim[0] << " " << RMDim[1] << " " << RMDim[2] << endl;

			ss << "X_COORDINATES " << RMDim[0] << " float" << endl;
			if( d == 0 )
			{
				ss << origin[0] + Rdel[0]*uber << endl;
			}
			else
			{
				for( long int i = 0; i < RMDim[0]; i++ )
					ss << origin[0] + Rdel[0]*i << " ";
				ss << endl;	
			}

			ss << "Y_COORDINATES " << RMDim[1] << " float" << endl;
			if( d == 1 )
			{
				ss << origin[1] + Rdel[1]*uber << endl;
			}
			else
			{
				for( long int i = 0; i < RMDim[1]; i++ )
					ss << origin[1] + Rdel[1]*i << " ";
				ss << endl;	
			}

			ss << "Z_COORDINATES " << RMDim[2] << " float" << endl;
			if( d == 2 )
			{
				ss << origin[2] + Rdel[2]*uber << endl;
			}
			else
			{
				for( long int i = 0; i < RMDim[2]; i++ )
					ss << origin[2] + Rdel[2]*i << " ";
				ss << endl;	
			}

			ss << "CELL_DATA " << (Rdims[id_1]-1)*(Rdims[id_2]-1) << endl;
			ss << "SCALARS isValid int" << endl;
			ss << "LOOKUP_TABLE default" << endl;

			long int count = 0;
			bool * dataP = data[d];

//Need to print in xyz order for vtk file.

			for( long int z_id = 0; z_id < RcMDim[2]; z_id++ )
			{	
				for( long int y_id = 0; y_id < RcMDim[1]; y_id++ )
				{	
					for( long int x_id = 0; x_id < RcMDim[0]; x_id++ )
					{

						if     ( d == 0 ) x_id = uber;
						else if( d == 1 ) y_id = uber;
						else if( d == 2 ) z_id = uber;
							
						long int id = x_id + (Rdims[0]-1) * ( y_id + (Rdims[1]-1) * z_id );  

						if( count%9 == 0 && count != 0 ){
							ss << endl;
						}
						ss << dataP[id] << " ";
						count++;

					}
				}
			}

			ss << endl;

			ss << "FIELD FieldData 1" << endl; 
			ss << "max_diff 1 " << (Rdims[id_1]-1)*(Rdims[id_2]-1) << " float" << endl;

			count = 0;
			double * diffP = diff[d];


			double max = 0;
			double min = (double)DBL_MAX;
			double avg = 0;
			long int counter = 0;

//Need to print in xyz order for vtk file.

			for( long int z_id = 0; z_id < RcMDim[2]; z_id++ )
			{	
				for( long int y_id = 0; y_id < RcMDim[1]; y_id++ )
				{	
					for( long int x_id = 0; x_id < RcMDim[0]; x_id++ )
					{

						if     ( d == 0 ) x_id = uber;
						else if( d == 1 ) y_id = uber;
						else if( d == 2 ) z_id = uber;
							
						long int id = x_id + (Rdims[0]-1) * ( y_id + (Rdims[1]-1) * z_id );  

						if( count%9 == 0 && count != 0 ){
							ss << endl;
						}

						double loc_diff = diffP[id];

						if( loc_diff > max )
						{
							max = loc_diff;
						}
						if( loc_diff != 0.0 ){
							if( loc_diff < min )
							{
								min = loc_diff;
							}
							avg += loc_diff;
							counter++;
						}
						ss << loc_diff << " ";
						count++;

					}
				}
			}

			if( counter == 0 )
			{
				avg = 0.0;
				min = 0.0;
			}
			else
			{
				avg = avg / counter;
				totAVG += avg;
				totSum++;
			}

/*
			cerr << "Dim:         " << ( (d == 0) ? 'X' : (d==1) ? 'Y' : 'Z' ) << endl;
			cerr << "Uber:        " << uber << endl;
			cerr << "NonZero Min: " << min << endl;
			cerr << "Max:         " << max << endl;
			cerr << "NonZero Avg: " << avg << endl;
			cerr << "Count:       " << counter << endl;
*/
			stringstream filename;
			filename << baseName << "_vtkOut_" << d << "_" << uber << ".vtk";
			string fname;
			filename >> fname;
			ofstream outFile( fname.c_str() );

			outFile << ss.rdbuf();

			outFile.close();
		}// uber
	}// d

	cerr << endl;
	cerr << "Total Avg: " << totAVG / ((totSum==0)?1:totSum) << endl;

}
