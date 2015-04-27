#include "Mesh.h"
#include "Particle.h"
#include "ParticleContainer.h"
#include "DomainMesh.h"

long int GetParticleCellID( Mesh* FineMesh, DomainMesh *UberMesh, Particle* part, double *MBB );
void ComputeFaceID( long int cellID, DomainMesh *UberMesh, Particle* part );
bool AcceptableFlow( long int cellID, DomainMesh *UberMesh, Mesh* FineMesh, Particle* part );
int AdvectParticleOnFlow( long int &cellID, Mesh *FineMesh, DomainMesh *UberMesh, Particle* part, double* MBB );
void AdvectParticleList( Mesh* FineMesh, DomainMesh* UberMesh, ParticleContainer* advectList, double endtime, double* MBB );
