#include "PrintVTK.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstring>

using namespace std;

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

	for( int d = 0; d < 3; d++ )
	{
		//Rdims from Flow_Mesh.C	
		//RMDim with 1 for ubers	
		//Rdel from Flow_Mesh.C	
		//Rpos from Flow_Mesh.C	
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

			stringstream filename;
			filename << baseName << "_vtkOut_" << d << "_" << uber << ".vtk";
			string fname;
			filename >> fname;
			ofstream outFile( fname.c_str() );

			outFile << ss.rdbuf();

			outFile.close();
		}// uber
	}// d
}
