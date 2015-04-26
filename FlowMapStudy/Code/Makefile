#CPP_COMPILER=g++
CPP_COMPILER=icpc

#OPLVL=-g
#OPLVL=-O0
#OPLVL=-O1
OPLVL=-O2
#OPLVL=-O3

COMPILE_FLAGS=-fopenmp

INCLUDE_DIRS=-IDomainMesh/ -IFlow/ -IMesh/ -IParticle/ -IParticleContainer/ -IPrintVTK/ -IProblem/

tokC=Problem/PA_Tokamak.C

domC=DomainMesh/DomainMesh.C
domH=DomainMesh/DomainMesh.h
domO=DomainMesh/DomainMesh.o

meshC=Mesh/Mesh.C
meshH=Mesh/Mesh.h
meshO=Mesh/Mesh.o

parC=Particle/Particle.C
parH=Particle/Particle.h
parO=Particle/Particle.o

pcnC=ParticleContainer/ParticleContainer.C
pcnH=ParticleContainer/ParticleContainer.h
pcnO=ParticleContainer/ParticleContainer.o

ptkC=PrintVTK/PrintVTK.C
ptkH=PrintVTK/PrintVTK.h
ptkO=PrintVTK/PrintVTK.o

floC=Flow/Flow.C
floH=Flow/Flow.h
floO=Flow/Flow.o

pntC=Flow/Point.C
pntH=Flow/Point.h
pntO=Flow/Point.o

avpC=Problem/AdvectParticles.C
avpH=Problem/AdvectParticles.h
avpO=Problem/AdvectParticles.o

OBJ=\
	$(ptkO) \
	$(avpO) \
	$(domO) \
	$(meshO) \
	$(pcnO) \
	$(floO) \
	$(parO) \
	$(pntO) 

EPSILON=0.0000001

STEPSIZE=0.001

VARS=-Depsilon=$(EPSILON) -DSTEPSIZE=$(STEPSIZE)

all: Tokamak 

Tokamak: $(tokC) $(OBJ)
	$(CPP_COMPILER) $(COMPILE_FLAGS) $(OPLVL) $(INCLUDE_DIRS) $(VARS) $(tokC) $(OBJ) -o bin/runTokamak.out

$(ptkO): $(ptkH) $(ptkC)
	$(CPP_COMPILER) $(COMPILE_FLAGS) $(OPLVL) $(INCLUDE_DIRS) -c $(ptkC) -o $(ptkO)

$(avpO): $(avpH) $(avpC) $(domO)
	$(CPP_COMPILER) $(COMPILE_FLAGS) $(OPLVL) $(INCLUDE_DIRS) $(VARS) -c $(avpC) $(domO) -o $(avpO)

$(domO): $(domH) $(domC) $(pcnO)
	$(CPP_COMPILER) $(COMPILE_FLAGS) $(OPLVL) $(INCLUDE_DIRS) $(VARS) -c $(domC) $(pcnO) -o $(domO)

$(pcnO): $(pcnH) $(pcnC) $(meshO)
	$(CPP_COMPILER) $(COMPILE_FLAGS) $(OPLVL) $(INCLUDE_DIRS) $(VARS) -c $(pcnC) $(meshO) -o $(pcnO) 

$(meshO): $(meshH) $(meshC) $(floO)
	$(CPP_COMPILER) $(COMPILE_FLAGS) $(OPLVL) $(INCLUDE_DIRS) $(VARS) -c $(meshC) $(floO) -o $(meshO)

$(floO): $(floH) $(floC) $(parO) 
	$(CPP_COMPILER) $(COMPILE_FLAGS) $(OPLVL) $(INCLUDE_DIRS) $(VARS) -c $(floC) $(parO) -o $(floO) 

$(parO): $(parH) $(parC) $(pntO) 
	$(CPP_COMPILER) $(COMPILE_FLAGS) $(OPLVL) $(INCLUDE_DIRS) $(VARS) -c $(parC) $(pntO) -o $(parO)

$(pntO): $(pntH) $(pntC)
	$(CPP_COMPILER) $(COMPILE_FLAGS) $(OPLVL) $(INCLUDE_DIRS) -c $(pntC) -o $(pntO)
	
clean:
	@rm */*.o bin/*.out

cvtk:
	@rm *.vtk

clear:
	@clear
	@clear
