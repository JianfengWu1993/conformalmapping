#include <fstream>
#include <vector>

#include "OBJFileReader.h"
#include "Solid.h"
#include "iterators.h"
#include "CTrait.h"
#include <math.h>
#include "Point.h"

using namespace MeshLib;

#define DELTA_T 1e-2
#define DELTA_E 1e-5


void starmap(Solid *cmesh)
{
	int count = 0;
	Point sum(0, 0, 0);
	for (SolidVertexIterator viter(cmesh); !viter.end(); ++viter)
	{
		Solid::tVertex vertex = *viter;
		sum += vertex->point();
		count++;		
	}
	Point center = sum / count;
	for (SolidVertexIterator viter(cmesh); !viter.end(); ++viter)
	{
		Solid::tVertex vertex = *viter;
		vertex->point() -= center;
		vertex->point() /= vertex->point().norm();
		v_n(vertex) = vertex->point();
	}
}

void gaussmap(Solid *cmesh)
{
	for (SolidVertexIterator viter(cmesh); !viter.end(); ++viter)
	{
		Solid::tVertex vertex = *viter;
		int count = 0;
		Point sum(0, 0, 0);
		for (VertexFaceIterator vfiter(vertex); !vfiter.end(); ++vfiter)
		{
			Solid::tFace face = *vfiter;
			sum += f_n(face);
			count++;
		}
		sum /= sum.norm();
		v_n(vertex) = sum;
		vertex->point() = v_n(vertex);
	}
}
double energy(Solid *cmesh) 
{
	double energy = 0;
	for (SolidEdgeIterator eiter(cmesh); !eiter.end(); ++eiter)
	{
		Solid::tEdge edge = *eiter;
		Solid::tVertex v1 = cmesh->edgeVertex1(edge);
		Solid::tVertex v2 = cmesh->edgeVertex2(edge);
		Point uv = v1->point() - v2->point();
		energy += uv.norm2();
	}
	return energy;
}

double energy_harmonic(Solid *cmesh)
{
	double energy = 0;
	for (SolidEdgeIterator eiter(cmesh); !eiter.end(); ++eiter)
	{
		Solid::tEdge edge = *eiter;
		Solid::tVertex v1 = cmesh->edgeVertex1(edge);
		Solid::tVertex v2 = cmesh->edgeVertex2(edge);
		Point uv = v1->point() - v2->point();
		energy += e_k(edge)*uv.norm2();
	}
	return energy;
}

void absderivative(Solid *cmesh)
{
	for (SolidVertexIterator viter(cmesh); !viter.end(); ++viter)
	{
		Solid::tVertex sv = *viter;
		Point gradient = Point(0, 0, 0);
		for (VertexVertexIterator vviter(sv); !vviter.end(); ++vviter)
		{
			Vertex * tv = *vviter;
			gradient += sv->point() - tv->point();
		}
		v_dv(sv) = gradient;
		v_abdv(sv) = v_dv(sv) - v_n(sv) * (v_dv(sv)*v_n(sv));
	}
}

void absderivative_harmonic(Solid *cmesh)
{
	for (SolidVertexIterator viter(cmesh); !viter.end(); ++viter)
	{
		Solid::tVertex sv = *viter;
		Point gradient = Point(0, 0, 0);
		for (VertexVertexIterator vviter(sv); !vviter.end(); ++vviter)
		{
			Vertex * tv = *vviter;
			double kuv = e_k(cmesh->vertexEdge(sv, tv));
			gradient += (sv->point() - tv->point())*kuv;
		}
		v_dv(sv) = gradient;
		v_abdv(sv) = v_dv(sv) - v_n(sv) * (v_dv(sv)*v_n(sv));
	}
}

void conjugategradient(Solid *cmesh)
{

	for (SolidVertexIterator viter(cmesh); !viter.end(); ++viter)
	{
		Solid::tVertex vertex = *viter;
		v_beta(vertex) = v_abdv(vertex)*v_abdv(vertex) / (v_fabdv(vertex)*v_fabdv(vertex));
		v_s(vertex) = v_s(vertex)*v_beta(vertex) - v_abdv(vertex);
	}
	double alpha = 1e-6;
	double energy = energy_harmonic(cmesh);
	for (SolidVertexIterator viter(cmesh); !viter.end(); ++viter)
	{
		Solid::tVertex v = *viter;
		v_mp(v) = v->point();
	}
	while (alpha<1e-2)
	{
		for (SolidVertexIterator viter(cmesh); !viter.end(); ++viter)
		{
			Solid::tVertex v = *viter;
			v->point() = v->point() + v_s(v)*alpha;
			v->point() /= v->point().norm();
		}
		double nenergy = energy_harmonic(cmesh);
		if (nenergy < energy)
		{
			for (SolidVertexIterator viter(cmesh); !viter.end(); ++viter)
			{
				Solid::tVertex v = *viter;
				v_mp(v) = v->point();
			}
		}
		alpha *= 2;
	}
	for (SolidVertexIterator viter(cmesh); !viter.end(); ++viter)
	{
		Solid::tVertex v = *viter;
		v->point() = v_mp(v);
	}



}
void updatemesh(Solid *cmesh)
{
	for (SolidVertexIterator viter(cmesh); !viter.end(); ++viter)
	{
		Solid::tVertex v = *viter;
		v->point() = v->point()- v_abdv(v)*DELTA_T;
		v->point() /= v->point().norm();
	}

}

void updatekuv(Solid * cmesh)
{
	double sum = 0;
	for (SolidEdgeIterator eiter(cmesh); !eiter.end(); ++eiter)
	{
		Solid::tEdge edge = *eiter;
		Point p1 = cmesh->edgeVertex1(edge)->point();
		Point p2 = cmesh->edgeVertex2(edge)->point();
		Point p3 = edge->halfedge(0)->he_next()->target()->point();

		double a, b;
		a = (p1 - p3)*(p2 - p3) / ((p1 - p3) ^ (p2 - p3)).norm() / 2;
		p3 = edge->halfedge(0)->he_sym()->he_next()->target()->point();
		b = (p1 - p3)*(p2 - p3) / ((p1 - p3) ^ (p2 - p3)).norm() / 2;
		e_k(edge) = a + b;
		sum += e_k(edge);
	}

}

void masscenter(Solid * cmesh)
{
	for (SolidFaceIterator fiter(cmesh); !fiter.end(); ++fiter)
	{
		Solid::tFace face = *fiter;

		Point p[3];
		int   i = 0;
		for (FaceVertexIterator fviter(face); !fviter.end(); ++fviter)
		{
			Vertex * v = *fviter;
			p[i++] = v->point();
		}

		Point n = (p[1] - p[0]) ^ (p[2] - p[0]);
		f_a(face) = n.norm() / 2.0;
	}
	for (SolidVertexIterator viter(cmesh); !viter.end(); ++viter) {
		Vertex * vertex = *viter;
		double area = 0;
		for (VertexFaceIterator vfiter(vertex); !vfiter.end(); ++vfiter) {
			Face* vf = *vfiter;
			area += f_a(vf);
		}
		v_a(vertex) = area / 3.0;
	}

	Point center = { 0, 0, 0 };
	double mass = 0;
	for (SolidVertexIterator viter(cmesh); !viter.end(); ++viter) {
		Vertex * vertex = *viter;
		center += vertex->point() * v_a(vertex);
		mass += v_a(vertex);
	}
	center /= mass;

	for (SolidVertexIterator viter(cmesh); !viter.end(); ++viter) {
		Vertex * vertex = *viter;
		vertex->point() -= center;
		vertex->point() /= vertex->point().norm();
	}

}

void main(int argc, char *argv[])
{
	// Read in the obj file
	Solid mesh;
	OBJFileReader of;
	std::ifstream in(argv[1]);
	of.readToSolid(&mesh, in);

	/******************* Put you subdivision processing here *********************/
	Solid nmesh;
	nmesh.add(mesh);
	//set trait
	for (SolidVertexIterator viter(&nmesh); !viter.end(); ++viter)
	{
		Solid::tVertex vertex = *viter;
		vertex->trait() = new CVertexTrait;
	}
	for (SolidFaceIterator fiter(&nmesh); !fiter.end(); ++fiter)
	{
		Solid::tFace face = *fiter;
		face->trait() = new CFaceTrait;
		Point fp[3]; 
		int i = 0;
		for (FaceVertexIterator fviter(face); !fviter.end(); ++fviter) {
			Vertex* fv = *fviter;
			fp[i++] = fv->point();
		}
		Point normal = (fp[1] - fp[0]) ^ (fp[2] - fp[0]);
		f_n(face) = normal / normal.norm();
	}
	for (SolidEdgeIterator eiter(&nmesh); !eiter.end(); ++eiter)
	{
		Solid::tEdge edge = *eiter;
		edge->trait() = new CEdgeTrait;
		e_k(edge) = 1.0;
	}
	//set edge kuv
	Solid * pnmesh = &nmesh;
	updatekuv(pnmesh);

	//Tuette mapping 
	starmap(pnmesh);

	//Tuette mapping step 2 to 5
	double te = energy(pnmesh);
	double oe = 1000;
	std::cout << "Tuette mapping start" << std::endl;
	std::cout << "Delta of Tuette energy and Tuette energy" << std::endl;
	while (oe-te > DELTA_E)
	{
		absderivative(pnmesh);
		updatemesh(pnmesh);
		oe = te;
		te = energy(pnmesh);
		std::cout << oe-te <<" " <<te<<std::endl;
	}


	//harmonic mapping with steepest descent
	//double he=energy_harmonic(pnmesh);

	//while (oe-he>DELTA_E)
	//{
	//	absderivative_harmonic(pnmesh);
	//	updatemesh(pnmesh);
	//	masscenter(pnmesh);
	//	oe = he;
	//	he = energy_harmonic(pnmesh);
	//	std::cout << oe - he <<" "<< he << std::endl;
	//}


	//harmonic mapping with conjugate gradient.
	std::cout << "Harmonic mapping start" << std::endl;
	std::cout << "Delta of Harmonic energy and Harmonic energy" << std::endl;
	absderivative_harmonic(pnmesh);
	updatemesh(pnmesh);
	for (SolidVertexIterator viter(pnmesh); !viter.end(); ++viter)
	{
		Solid::tVertex vertex = *viter;
		v_fabdv(vertex) = v_abdv(vertex);
		v_s(vertex) = v_abdv(vertex);
	}
	double he = energy_harmonic(pnmesh);
	while (oe - he > DELTA_E)
	{
		absderivative_harmonic(pnmesh);
		conjugategradient(pnmesh);
		oe = he;
		he = energy_harmonic(pnmesh);
		std::cout << oe - he << " " << he << std::endl;
	}

	// Write out the resultant obj file
	int vObjID = 1;
	std::map<int, int> vidToObjID;

	std::ofstream os(argv[2]);

	SolidVertexIterator iter(&nmesh);

	for(; !iter.end(); ++iter)
	{
		Vertex *v = *iter;
		Point p = v->point();
		os << "v " << p[0] << " " << p[1] << " " << p[2] << std::endl;
		vidToObjID[v->id()] = vObjID++;
	}
	os << "# " << (unsigned int)nmesh.numVertices() << " vertices" << std::endl;

	float u = 0.0, v = 0.0;
	for(iter.reset(); !iter.end(); ++iter)
	{
		Vertex *vv = *iter;
		std::string key( "uv" );
		std::string s = Trait::getTraitValue (vv->string(), key );
		if( s.length() > 0 )
		{
			sscanf( s.c_str (), "%f %f", &u, &v );
		}
		os << "vt " << u << " " << v << std::endl;
	}
	os << "# " << (unsigned int)nmesh.numVertices() << " texture coordinates" << std::endl;

	SolidFaceIterator fiter(&nmesh);
	for(; !fiter.end(); ++fiter)
	{
		Face *f = *fiter;
		FaceVertexIterator viter(f);
		os << "f " ;
		for(; !viter.end(); ++viter)
		{
			Vertex *v = *viter;
			os << vidToObjID[v->id()] << "/" << vidToObjID[v->id()] << " ";
		}
		os << std::endl;
	}
	os.close();
	std::cout << "Press Enter to exit." << std::endl;
	getchar();
}