
#include <assert.h>
#include "iterators.h"

#include "SolidDelegate.h"
#include "EditTrait.h"

using namespace MeshLib;

Vertex * SolidDelegate::edgeSplit(Solid *pS, Edge *e)
{
	return pS->edgeSplit(e);
}

bool SolidDelegate::collapsable (Solid *pS, Edge * e)
{
	return pS->collapsable (e->halfedge (0));
}

void SolidDelegate::collapseEdge (Solid *pS, Edge * e)
{
	pS->collapseEdge (e);
}
Vertex * SolidDelegate::createVertex( Solid *pS, int id )
{
	return pS->createVertex(id);
}

Face * SolidDelegate::createFace( Solid * pS, int * v, int id)
{
	return pS->createFace(v,id);
}

void SolidDelegate::removeEdge( Solid * pS, Edge * e)
{
	Vertex * v1 = pS->idVertex (e->ekey ().s ());
	Vertex * v2 = pS->idVertex (e->ekey ().t ());
	if(v1!= NULL && v1->halfedge () == NULL)
	{
		removeVertex(pS,v1);
		delete v1;
	}
	if(v2 != NULL && v2->halfedge () == NULL)
	{
		removeVertex(pS,v2);
		delete v2;
	}
	pS->m_edges.remove(e);
}

void SolidDelegate::removeFace (Solid * pS, Face * f) 
{
	int i;
	HalfEdge * he = f->halfedge ();
	for(i=0; i<3; i++)
	{
		he = he->he_next ();
		Vertex * v = he->target ();
		if(v->halfedge () == he)
		{
			for(SolidEdgeIterator eiter(pS); !eiter.end (); ++eiter)
			{
				Edge * e = *eiter;
				HalfEdge * hhe = e->halfedge (0);
				if(hhe->target () == v && hhe != he)
				{
					v->halfedge () = hhe;
					break;
				}
				else if (e->halfedge (1) != NULL && e->halfedge (1)->target () == v && e->halfedge (1) != he)
				{
					v->halfedge () = e->halfedge (1);
					break;
				}
			}
		}
		if(v->halfedge () == he)
			v->halfedge () = NULL;
	}
	for(i=0; i<3; i++)
	{
		HalfEdge * nhe = he->he_next ();
	//	Vertex * v = he->target ();
		Edge * e = he->edge ();
		if(e->halfedge (1) == NULL)
		{
			removeEdge(pS, e);
			delete e;
			delete he;
		}
		else
		{
			if(e->halfedge (0) == he)
				e->halfedge (0) = e->halfedge (1);
			delete he;
			e->halfedge(1) = NULL;
		}
		he = nhe;
	}

	pS->m_faces.remove(f);
}

void SolidDelegate::removeVertex( Solid *pS, Vertex * v)
{
	pS->m_verts.remove(v);
}



//
void SolidDelegate::removeHalfedge(Solid *pS,  HalfEdge * he )
{
	//he->edge
	Edge * edge = he->edge();

	if( edge->halfedge(0) == he )
	{
		edge->halfedge(0) = edge->halfedge(1);
		edge->halfedge(1) = NULL;
	}
	else
	{
		if( edge->halfedge(1) != he ) throw TopologyException("SolidDelegate");
		edge->halfedge(1) = NULL;
	}

	
	//he->face
	if( he->face()->halfedge() == he )
	{
		he->face()->halfedge() = NULL;
	}

	if( he->trait() != NULL )
	{
		Trait::clear(  he->trait() );
	}
	//
	if( edge->halfedge(0) == NULL && edge->halfedge(1) == NULL )
	{
		if( edge->trait() != NULL )
		{
			Trait::clear( edge->trait() );
		}
		pS->m_edges.remove( edge );
		delete edge;
		edge = NULL;
	}
	delete he;
};

void SolidDelegate::removeFace2(Solid *pS,  Face * f )
{
	if( f == NULL ) return;

	HalfEdge * helist[3];

	int i =0;
	for( FaceHalfedgeIterator fhiter(f); !fhiter.end(); ++fhiter)
	{
		HalfEdge * he = *fhiter;
		helist[i] = he;
		i++;
	}
	
	for( i =0; i < 3; i++)
	{
		Vertex *vt = helist[i]->target();
		if( vt->halfedge() == helist[i])
		{
			for( ListIterator<HalfEdge> liter( *v_edit_halfedge(vt)); !liter.end(); ++liter)
			{
				HalfEdge *he = *liter;
				if( he != helist[i])
				{
					vt->halfedge() = he;
					break;
				}
			}
		}
	}
	for( i = 0; i < 3; i ++ )
	{
		removeHalfedge(pS, helist[i] );
	}


	pS->m_faces.remove( f );

	if( f->trait() != NULL )
	{
		Trait::clear( f->trait() );
	}
	delete f;
}



void SolidDelegate::destroyEdge (Solid *pS, Edge * e){
	pS->destroyEdge (e);
}

void SolidDelegate::destroyVertex (Solid *pS, Vertex * v)
{
	pS->destroyVertex (v);
}

void SolidDelegate::destroyFace (Solid *pS, Face * f)
{
	pS->destroyFace (f);
}

void SolidDelegate::destroyHalfEdge (Solid * pS, HalfEdge * he)
{
	pS->destroyHalfEdge (he);
}

void SolidDelegate::collapseEdgeVertex( Solid *pS, Edge *edge, Vertex * vertex )
{
	pS->collapseEdgeVertex( edge, vertex );
}

