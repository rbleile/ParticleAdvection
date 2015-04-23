all: euler

icc:
	icc PA_Flow.C -g -Wall -o patedr

euler: FlowMap.h Mesh.h
	g++ PA_Flow.C -o paFlowEuler
