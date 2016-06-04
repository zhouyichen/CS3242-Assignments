#include "Mesh.h"
#include "Eigen/Dense"
#include <map>
#include <unordered_map>

using namespace std;
using namespace Eigen;

Mesh::Mesh(const char* filename){
	loadMF(filename);
}

map<int, Point> allVertices; // to store all the active vertices
map<int, Matrix4d> ErrorMatrix; // to store all the error matrix for each vertex
map<int, Triangle> faces; // to store all the active faces/triangles

unordered_map<int, pair<double, Vector3d>> edgesErrors; //Used as a hashmap to find active edges
map<double, pair<int, int>> errorsEdges; //Used as a heap to get the edge to be collapsed

void Mesh::loadMF(const char* filename){
	std::ifstream infile;
	infile.open(filename, std::ios::in);
	std::string strbuff;
	int index = 1;
	while(std::getline(infile,strbuff)){
		std::stringstream ss;
		ss<<strbuff;
		char type;
		ss>>type;
		if(type=='v'){
			Vertex v;
			ss>>v.x>>v.y>>v.z;
			new Point(v.x, v.y, v.z, index);
			index++;
		}
		else if(type=='f'){
			Triangle t;
			ss>>t.p1>>t.p2>>t.p3;
			t.id = index;
			t.setTriangle();
			index++;
		}
	}
	infile.close();
}

void Point::setPoint (double x, double y, double z, int i){
	id = i;
	v[0] = x;
	v[1] = y;
	v[2] = z;
	ErrorMatrix.emplace(id, Matrix4d::Zero());
	allVertices.emplace(id, *this);
}
	
void Point:: setVector(Vector3d new_v){
	this ->v = new_v;
}

void Triangle:: setTriangle(){
	Point *p_1 = &allVertices[p1];
	Point *p_2 = &allVertices[p2];
	Point *p_3 = &allVertices[p3];
	
	// Add faces index to each vertex
	// Add edges index to each vertex;
	p_1->faces.emplace(id);
	p_1->neighbors.emplace(p2);
	p_1->neighbors.emplace(p3);
	p_2->faces.emplace(id);
	p_2->neighbors.emplace(p1);
	p_2->neighbors.emplace(p3);
	p_3->faces.emplace(id);
	p_3->neighbors.emplace(p2);
	p_3->neighbors.emplace(p1);
	
	// Calculate the plane equation and matrix
	Vector3d normal = ((p_1->v - p_2->v).cross(p_1->v - p_3->v)).normalized();
	double d = -1*normal.dot(p_1->v);
	plane << normal, d;
	mat = plane * (plane.transpose());
	
	// Add the matrix to each vertex in the ErrorMatrix
	ErrorMatrix[p1] += mat;
	ErrorMatrix[p2] += mat;
	ErrorMatrix[p3] += mat;
	
	// store the triangle in map
	faces.emplace(id, *this);
}

bool Triangle:: contains(int n){
	return p1==n || p2==n || p3==n;
}

// Get the minimum error for the egde between p1 and p2, with the vertex achieving that value
pair<double, Vector3d> getMinErrorNewVertex (int i1, int i2){
	
	Vector4d a(0, 0, 0, 1);
	Matrix4d errorMat = ErrorMatrix[i1] + ErrorMatrix[i2];
	
	Matrix4d Q = errorMat;
	Q.row(3) << 0, 0, 0, 1;
	
	Vector4d solution = Q.colPivHouseholderQr().solve(a);
	Vector3d newV(solution.x(), solution.y(), solution.z());
	
	double minError = solution.transpose() * errorMat * solution;
	return pair<double, Vector3d>(minError,newV);
}

void initialiseEdges(){
	// Iterate through all the triangles
	int a, b, c;
	int temp;
	for (auto t: faces){
		a = t.second.p1;
		b = t.second.p2;
		c = t.second.p3;
		if (b > c){
			temp = c; c = b; b = temp;
		}
		if (a > b){
			temp = b; b = a; a = temp;
			if (b > c){
				temp = c; c = b; b = temp;
			}
		}
		// Add 3 edges to the edgeErrors and errorEdges
		pair<double, Vector3d> ab = getMinErrorNewVertex(a,b);
		pair<double, Vector3d> bc = getMinErrorNewVertex(b,c);
		pair<double, Vector3d> ac = getMinErrorNewVertex(a,c);

		edgesErrors.emplace(a*1000000+b, ab);
		edgesErrors.emplace(b*1000000+c, bc);
		edgesErrors.emplace(a*1000000+c, ac);
		errorsEdges.emplace(ab.first, pair<int,int>(a,b));
		errorsEdges.emplace(bc.first, pair<int,int>(b,c));
		errorsEdges.emplace(ac.first, pair<int,int>(a,c));
	}
}

void Mesh::collapseEdge(){
	
	// get the edge with min error and the new position Vector3ds
	map<double, pair<int, int>>::iterator it = errorsEdges.begin();
	
	double error = it->first;
	pair<int, int> edge = it->second;
	errorsEdges.erase(it);
	
	int v1 = edge.first;
	int v2 = edge.second;

	auto edgeError = edgesErrors.find(v1*1000000+v2);
	if (edgeError == edgesErrors.end()){
		return;
	}
	else if (edgeError->second.first != error){
		return;
	}
	
	Point* p1 = &allVertices[v1];
	Point* p2 = &allVertices[v2];
	
	// Update v1 to the new position
	Vector3d newP = edgeError->second.second;
	p1 -> setVector(newP);
	
	// Update the matrix of v1 in ErrorMatrix M(v1) = M(v1)+M(v2)
	ErrorMatrix[v1] += ErrorMatrix[v2];
	
	// For all the faces incident with v2
	for (auto it : p2->faces){
		// if it incident with v1, remove it from v1,v2 and faces
		Triangle *t = &faces[it];
		if (t->contains(v1)){
			p1 ->faces.erase(it);
			// Remove the face from the third vertex forming triangle with v1 and v2
			if (t->p1 != v2 && t->p1 != v1){
				allVertices[t->p1].faces.erase(it);
			}
			else if (t->p2 != v2 && t->p2 != v1){
				allVertices[t->p2].faces.erase(it);
			}
			else {
				allVertices[t->p3].faces.erase(it);
			}
			faces.erase(it);
		}
		// else, replace v2 by v1, add the face to v1
		else{
			p1->faces.emplace(it);
			if (t->p1 == v2){t->p1 = v1;}
			else if (t->p2 == v2){t->p2 = v1;}
			else {t->p3 = v1;}
		}
	}
	// For all the neighbors of v2 (including v1)
	for (auto it: p2->neighbors){
		Point *p = &allVertices[it];
		// remove v2 from the neighbors of v
		p->neighbors.erase(v2);
		// remove the pair (v,v2) from the edgeErrors (including (v1, v2))
		if (it < v2){
			edgesErrors.erase(it*1000000+v2);
		}
		else{
			edgesErrors.erase(v2*1000000+it);
		}
		// if v is not v1, connect it with v1
		if (it != v1){
			p1->neighbors.emplace(it);
			p->neighbors.emplace(v1);
		}
	}
	// Update all the edges incident with v1
	for (auto it: p1 ->neighbors){
		pair<double, Vector3d> newError = getMinErrorNewVertex(it, v1);
		if (it < v1){
			edgesErrors[it*1000000+v1] = newError;
			errorsEdges[newError.first] = pair<int, int>(it, v1);
		}
		else{
			edgesErrors[v1*1000000+it] = newError;
			errorsEdges[newError.first] = pair<int, int>(v1, it);
		}
	}
	// Remove v2 from allVertices
	allVertices.erase(v2);
}


void Mesh::writeMF(const char* filename){
	std::ofstream outfile;
	outfile.open(filename, std::ios::out);
	std::string strbuff;
	
	map<int, int> id_index;
	int index = 1;
	for (auto it: allVertices){
		Vector3d v = it.second.v;
		id_index.emplace(it.second.id, index);
		outfile<<"v "<<v[0]<<" "<<v[1]<<" "<<v[2]<<std::endl;
		index ++;
	}
	for (auto it: faces){
		Triangle *t = &it.second;
		outfile<<"f "<<id_index[t->p1]<<" "<<id_index[t->p2]<<" "<<id_index[t->p3]<<std::endl;
	}
	outfile.close();
}

void Mesh::simplifyMesh(const char* input, const char* output, int faceCnt){
	// you may assume inputs are always valid
	loadMF(input);
	std::cout<<"Original face count: "<<faces.size()<<std::endl;
	std::cout<<"Original vertices count: "<< allVertices.size() <<endl;
	
	initialiseEdges();
	cout<<"Original edges count: "<< edgesErrors.size() <<endl;
	
	while (faces.size()>faceCnt){
		collapseEdge();
	}

	std::cout<<"Final vertices count: "<< allVertices.size() <<endl;
	cout<<"Final edges count: "<< edgesErrors.size() <<endl;
	cout<<"Final face count: "<< faces.size() <<endl;
	writeMF(output);
}
