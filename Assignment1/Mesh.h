#ifndef _MESH_H_
#define _MESH_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "Eigen/Dense"

#include <set>
#include <map>
using namespace std;
using namespace Eigen;

typedef unsigned int uint;

typedef struct {
	//3d coordinates
	float x;
	float y;
	float z;
} Vertex;

struct Point{
	Vector3d v;
	int id;
	void setPoint(double x, double y, double z, int i);
	Point(){};
	Point(double x, double y, double z, int i){setPoint(x,y,z, i);};
	void setVector(Vector3d new_v);
	set<int> neighbors;
	set<int> faces;
};

struct Triangle{
	int p1;
	int p2;
	int p3;
	int id;
	Vector4d plane;
	Matrix4d mat;
	void setTriangle();
	Triangle(){};
	bool contains(int n);
};

class Mesh{
public:
	Mesh() {};
	Mesh(const char*);
	//load a Mesh from .mesh file
	void loadMF(const char*);
	//write a Mesh to .mesh file (no header)
	void writeMF(const char*);
	//simplify a mesh
	void simplifyMesh(const char* input, const char* output, int faceCnt);
	void collapseEdge();
};
#endif