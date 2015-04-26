#ifndef POINT
#define POINT

#if 0
extern const double epsilon;
extern const int PointsToFaces[8][3];
extern const int FacesToPoints[6][4];
extern const int EdgesToFaces[12][2];
extern const int FacesToEdges[6][4];
extern const int FaceToFace[6];
extern const int BoundsToFace[6];
#endif

class Point
{
  public:
	double x;
	double y;
	double z;
	double t;

	Point();

	Point( double ix, double iy, double iz );

	Point( double ix, double iy, double iz, double it );

	void setPoint( double ix, double iy, double iz );

	void setPoint( double ix, double iy, double iz, double it );
		
};

#endif

