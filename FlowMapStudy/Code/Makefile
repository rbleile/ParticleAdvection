all: clear euler1 euler2

euler1: PA_Flow.C Flow_Mesh.o Mesh.o ParticleList.o Particle.o FlowMap.o PrintVTK.o
	g++ -O2 PA_Flow.C Flow_Mesh.o Mesh.o ParticleList.o Particle.o FlowMap.o PrintVTK.o -o paCustom.out

euler2: PA_Tokamak.C Flow_Mesh.o Mesh.o ParticleList.o Particle.o FlowMap.o PrintVTK.o
	g++ -O2 PA_Tokamak.C Flow_Mesh.o Mesh.o ParticleList.o Particle.o FlowMap.o PrintVTK.o -o paTokamak.out

PrintVTK.o: PrintVTK.h PrintVTK.C Flow_Mesh.h
	g++ -O2 -c PrintVTK.C -o PrintVTK.o

FlowMap.o: FlowMap.C FlowMap.h
	g++ -O2 -c FlowMap.C -o FlowMap.o

Particle.o: MeshParticle.C MeshParticle.h FlowMap.h FlowMap.o
	g++ -O2 -c MeshParticle.C -o Particle.o

ParticleList.o: MeshParticleList.C MeshParticleList.h MeshParticle.h Particle.o
	g++ -O2 -c MeshParticleList.C -o ParticleList.o

Mesh.o: Mesh.C Mesh.h MeshParticle.h MeshParticleList.h ParticleList.o Particle.o
	g++ -O2 -c Mesh.C -o Mesh.o

Flow_Mesh.o: Flow_Mesh.C Flow_Mesh.h FlowMap.h Mesh.h Mesh.o
	g++ -O2 -c Flow_Mesh.C -o Flow_Mesh.o 

debug: clear euler1d euler2d

euler1d: PA_Flow.C Flow_Meshd.o Meshd.o ParticleListd.o Particled.o FlowMapd.o PrintVTKd.o
	g++ -g PA_Flow.C Flow_Meshd.o Meshd.o ParticleListd.o Particled.o FlowMapd.o PrintVTKd.o -o paCustom.out

euler2d: PA_Tokamak.C Flow_Meshd.o Meshd.o ParticleListd.o Particled.o FlowMapd.o PrintVTKd.o
	g++ -g PA_Tokamak.C Flow_Meshd.o Meshd.o ParticleListd.o Particled.o FlowMapd.o PrintVTKd.o -o paTokamak.out

PrintVTKd.o: PrintVTK.h PrintVTK.C Flow_Mesh.h
	g++ -g -c PrintVTK.C -o PrintVTKd.o

FlowMapd.o: FlowMap.C FlowMap.h
	g++ -g -c FlowMap.C -o FlowMapd.o

Particled.o: MeshParticle.C MeshParticle.h FlowMap.h FlowMapd.o
	g++ -g -c MeshParticle.C -o Particled.o

ParticleListd.o: MeshParticleList.C MeshParticleList.h MeshParticle.h Particled.o
	g++ -g -c MeshParticleList.C -o ParticleListd.o

Meshd.o: Mesh.C Mesh.h MeshParticle.h MeshParticleList.h ParticleListd.o Particled.o
	g++ -g -c Mesh.C -o Meshd.o

Flow_Meshd.o: Flow_Mesh.C Flow_Mesh.h FlowMap.h Mesh.h Meshd.o
	g++ -g -c Flow_Mesh.C -o Flow_Meshd.o 

clean:
	@rm *.o *.out

cvtk:
	@rm *.vtk

clear:
	@clear
	@clear
