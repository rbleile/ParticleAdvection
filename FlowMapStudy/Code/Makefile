CPP_COMPILER=g++
#CPP_COMPILER=mpic++
#CPP_COMPILER=icpc

#OPLVL=-g
#OPLVL=-O0
#OPLVL=-O1
#OPLVL=-O2
OPLVL=-O3

#COMPILE_FLAGS=-openmp
#COMPILE_FLAGS= -Wall
COMPILE_FLAGS=-fopenmp

INCLUDE_DIRS=-IDomainMesh/ -IFlow/ -IMesh/ -IParticle/ -IParticleContainer/ -IPrintVTK/ -IProblem/

tokC=Problem/Tokamak.C
stokC=Problem/Tokamak_small.C

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

ftleC=Problem/FTLE.C
ftleH=Problem/FTLE.h
ftleO=Problem/FTLE.o

OBJ=\
	$(ptkO) \
	$(ftleO) \
	$(avpO) \
	$(domO) \
	$(meshO) \
	$(pcnO) \
	$(floO) \
	$(parO) \
	$(pntO) 

EPSILON=0.000001
DO_MPI=0

VARS=-Depsilon=$(EPSILON) -DDO_MPI=$(DO_MPI)

all: DIRS Tokamak TokamakS

DIRS:
	@mkdir bin -p

Tokamak: $(tokC) $(OBJ)
	$(CPP_COMPILER) $(COMPILE_FLAGS) $(OPLVL) $(INCLUDE_DIRS) $(VARS) $(tokC) $(OBJ) -o bin/runTokamak.out

TokamakS: $(stokC) $(OBJ)
	$(CPP_COMPILER) $(COMPILE_FLAGS) $(OPLVL) $(INCLUDE_DIRS) $(VARS) $(stokC) $(OBJ) -o bin/runTokamak_small.out

$(ptkO): $(ptkH) $(ptkC)
	$(CPP_COMPILER) $(COMPILE_FLAGS) $(OPLVL) $(INCLUDE_DIRS) $(VARS) -c $(ptkC) -o $(ptkO)

$(ftleO): $(ftleH) $(ftleC) $(parO)
	$(CPP_COMPILER) $(COMPILE_FLAGS) $(OPLVL) $(INCLUDE_DIRS) -c $(ftleC) -o $(ftleO)

$(avpO): $(avpH) $(avpC) $(domO) $(meshO) $(pcnO) $(parO)
	$(CPP_COMPILER) $(COMPILE_FLAGS) $(OPLVL) $(INCLUDE_DIRS) $(VARS) -c $(avpC) -o $(avpO)

$(domO): $(domH) $(domC) $(parO) $(floO) $(meshO)
	$(CPP_COMPILER) $(COMPILE_FLAGS) $(OPLVL) $(INCLUDE_DIRS) $(VARS) -c $(domC) -o $(domO)

$(pcnO): $(pcnH) $(pcnC) $(meshO) $(parO) $(floO)
	$(CPP_COMPILER) $(COMPILE_FLAGS) $(OPLVL) $(INCLUDE_DIRS) $(VARS) -c $(pcnC) -o $(pcnO) 

$(meshO): $(meshH) $(meshC) $(parO) $(floO)
	$(CPP_COMPILER) $(COMPILE_FLAGS) $(OPLVL) $(INCLUDE_DIRS) $(VARS) -c $(meshC) -o $(meshO)

$(floO): $(floH) $(floC) $(pntO)
	$(CPP_COMPILER) $(COMPILE_FLAGS) $(OPLVL) $(INCLUDE_DIRS) $(VARS) -c $(floC) -o $(floO) 

$(parO): $(parH) $(parC) $(pntO) 
	$(CPP_COMPILER) $(COMPILE_FLAGS) $(OPLVL) $(INCLUDE_DIRS) $(VARS) -c $(parC) -o $(parO)

$(pntO): $(pntH) $(pntC)
	$(CPP_COMPILER) $(COMPILE_FLAGS) $(OPLVL) $(INCLUDE_DIRS) -c $(pntC) -o $(pntO)
	
clean:
	@rm */*.o bin/*.out

cvtk:
	@rm *.vtk bin/*.vtk

clear:
	@clear
	@clear

count:
	wc -l `find *` || true
