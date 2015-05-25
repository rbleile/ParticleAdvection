#include <Particle.h>

double dist( Particle p1, Particle p2 );
double compute_FTLE(	int x_id, int y_id, int z_id, 
						const long int nx, const long int ny, const long int nz, 
						Particle* inputs, Particle* outputs, double endTime );
