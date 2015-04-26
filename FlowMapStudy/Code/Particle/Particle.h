#include "Point.h"
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

	double getStepSize();
	void setStepSize( double ss);
	void setFID( double* bbox);

};

#endif

