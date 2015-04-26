#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstring>

#ifndef PRINT_VTK_H
#define PRINT_VTK_H

using std::string;

class VTKFilePrinter
{
  private:
	int Udims[3];
	int Fdims[3];
	double Udel[3];
	double Fdel[3];
	double origin[3];
	bool **data;
	string baseName;
	
  public:
	
	VTKFilePrinter( int u[], int f[], double ud[], double fd[], double x[], bool ** d, string s )
	{
		Udims[0] = u[0];
		Udims[1] = u[1];
		Udims[2] = u[2];
		Fdims[0] = f[0];
		Fdims[1] = f[1];
		Fdims[2] = f[2];

		Udel[0] = ud[0];
		Udel[1] = ud[1];
		Udel[2] = ud[2];
		Fdel[0] = fd[0];
		Fdel[1] = fd[1];
		Fdel[2] = fd[2];

		origin[0] = x[0];
		origin[1] = x[1];
		origin[2] = x[2];

		data = d;

		baseName = s;
	}

	void printVtkDs();

};

#endif

