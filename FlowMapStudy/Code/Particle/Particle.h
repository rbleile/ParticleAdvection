#include <Point.h>
#include <cmath>
#include <ostream>
#include <iostream>
using std::ostream;
using std::cerr;
using std::endl;
using std::fabs;

/*
extern const double epsilon;
extern const int PointsToFaces[8][3];
extern const int FacesToPoints[6][4];
extern const int EdgesToFaces[12][2];
extern const int FacesToEdges[6][4];
extern const int FaceToFace[6];
extern const int BoundsToFace[6];
*/

#ifndef MESH_PARTICLE
#define MESH_PARTICLE

class Particle : public Point
{
  private:
	double stepSize;
  public:
	int FID;

	Particle(){};
	Particle( Particle *p )
	{
		x = p->x;
		y = p->y;
		z = p->z;
		t = p->t;
		stepSize = p->getStepSize();
		FID = p->FID;
	}
	double getStepSize();
	void setStepSize( double ss);
	void setFID( double* bbox);

};

#endif

