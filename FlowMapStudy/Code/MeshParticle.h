#include "FlowMap.h"

#ifndef MESH_PARTICLE
#define MESH_PARTICLE

class Particle : public FMPoint
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

