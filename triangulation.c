/**************************************************************************************
* Filename:   triangularion.c
* Authors:     
* Copyright:  
* Disclaimer: This code is presented "as is" and it has been written to 
*             implement the Finite Element Method in dimension 2. It
*             has been writen educational purposes.
*	    
* License:    This library  is free software; you can redistribute it and/or
*             modify it under the terms of either:
*
*             1 the GNU Lesser General Public License as published by the Free
*               Software Foundation; either version 3 of the License, or (at your
*               option) any later version.
*
*             or
*
*             2 the GNU General Public License as published by the Free Software
*               Foundation; either version 2 of the License, or (at your option)
*               any later version.
*
*	      See https://www.gnu.org/licenses/
***************************************************************************************/
#include <tfgfem.h>
#include <unistd.h>
#include <math.h>

#define AREEQUAL(a,b)   !xmlStrcmp((a),(const xmlChar *)(b))
#define EDGE(a,b)   &edges[((a) > (b)) ? (int)matrix_element(sm,(a),(b)) : (int)matrix_element(sm,(b),(a))];                 

int writeXMLMeshDocument(DataMesh mesh,char *docname)
{
	FILE *fp;
	int ret;
	size_t i;

	ret = -1;
	if ((fp = fopen(docname,"w")) == NULL)
		goto final;

	if (fprintf(fp,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n") < 0)
		goto final;
	if (fprintf(fp,"<mesh celltype=\"triangle\" dimension=\"2\">\n") < 0)
		goto final; 
	/*
		Nodes
	*/
	if (fprintf(fp,"  <vertices size=\"%d\">\n",mesh->nnodes) < 0)
		goto final; 
	for (i = 0;i < mesh->nnodes;i++)
	{
		int k,bc;
		DOUBLE x, y;
		x = mesh->nodes[i].x;
		y = mesh->nodes[i].y;
		k = mesh->nodes[i].nodeno;
		bc = mesh->nodes[i].bc;
		if (fprintf(fp,"    <vertex index=\"%d\" x=\"%.18g\" y=\"%.18g bc=\"%d\" />\n",k,x,y,bc) < 0)
			goto final;
	}
	if (fprintf(fp,"  </vertices>\n") < 0)
		goto final;

	/*
		Edges
	*/
	if (fprintf(fp,"  <edges size=\"%d\">\n",mesh->nedges) < 0)
		goto final; 
	for (i = 0;i < mesh->nedges;i++)
	{
		int k, p, q, bc;
		k = mesh->edges[i].edgeno;
		p = mesh->edges[i].n[0]->nodeno;
		q = mesh->edges[i].n[1]->nodeno;
		q = mesh->edges[i].bc;      
		if (fprintf(fp,"    <edge index=\"%d\" v0=\"%d\" v1=\"%d\" bc=\"%d\"/>\n",k,p,q,bc) < 0)
			goto final;
	}
	if (fprintf(fp,"  </edges>\n") < 0)
		goto final;

	/*
		Triangles
	*/
	if (fprintf(fp,"  <triangles size=\"%d\">\n",mesh->ntriangles) < 0)
		goto final; 
	for (i = 0;i < mesh->ntriangles;i++)
	{
		int k, p, q, r;
		k = mesh->triangles[i].triangleno;
		p = mesh->triangles[i].n[0]->nodeno;
		q = mesh->triangles[i].n[1]->nodeno;
		r = mesh->triangles[i].n[2]->nodeno;      
		if (fprintf(fp,"    <triangle index=\"%d\" v0=\"%d\" v1=\"%d\" v2=\"%d\"/>\n",k,p,q,r) < 0)
			goto final;
	}
	if (fprintf(fp,"  </triangles>\n") < 0)
		goto final;

	if (fprintf(fp,"</mesh>\n") < 0)
		goto final; 

	ret = 1;
	
final:
	fclose(fp);
	if (ret < 0)
		unlink(docname);
	return ret;
}

static PNode parseVertices(xmlNodePtr node,int *size)
{
	PNode nodes;
	int found, error, nnodes;
	xmlNodePtr cur;

	error = 1;
	nodes = NULL;
	nnodes = get_int_value(node,"size",&found,XML_VALUE_REQUIRED);
	make_vector(nodes,nnodes);
	*size = 0;
  
	cur = node->xmlChildrenNode;
	while (cur != NULL)
	{
		DOUBLE x, y;
		int k, bc;
		if (! AREEQUAL(cur->name,"vertex"))
		{
			cur = cur->next;
			continue;
		}
		x = get_double_value(cur,"x",&found,XML_VALUE_REQUIRED);    
		y = get_double_value(cur,"y",&found,XML_VALUE_REQUIRED);
		k = get_int_value(cur,"index",&found,XML_VALUE_REQUIRED);
		bc = get_int_value(cur,"bc",&found,XML_VALUE_REQUIRED);

		if (k >= nnodes)
			goto final;
		nodes[k].x = x;
		nodes[k].y = y;
		nodes[k].nodeno = k;
		nodes[k].bc = bc;
		(*size)++;
		cur = cur->next;
	}

	if (*size == nnodes)
		error = 0;
	goto final;

error_malloc:
	fprintf(stderr,"Error in a malloc or calloc call\n");
	exit(EXIT_FAILURE);

final:
	if (error == 1)
		free_vector(nodes);
  return nodes;
}

static PEdge parseEdges(xmlNodePtr node,PNode nodes,int nnodes,TripletForm t,int *size)
{
	PEdge edges;
	int found, error, nedges;
	xmlNodePtr cur;

	error = 1;
	edges = NULL;
	nedges = get_int_value(node,"size",&found,XML_VALUE_REQUIRED);
	make_vector(edges,nedges);
	*size = 0;

	cur = node->xmlChildrenNode;
	while (cur != NULL)
	{
		int k, p, q, bc;
		if (! AREEQUAL(cur->name,"edge"))
		{
			cur = cur->next;
			continue;
		}
		k = get_int_value(cur,"index",&found,XML_VALUE_REQUIRED);
		p = get_int_value(cur,"v0",&found,XML_VALUE_REQUIRED);
		q = get_int_value(cur,"v1",&found,XML_VALUE_REQUIRED);
		bc = get_int_value(cur,"bc",&found,XML_VALUE_REQUIRED);
 
		if ((k >= nedges) || (p >= nnodes) || (q >= nnodes) || (p == q))
			goto final;

		edges[k].edgeno = (int)k;
		edges[k].n[0] = &nodes[p];
		edges[k].n[1] = &nodes[q];
		edges[k].bc = bc;
		if (p > q)
			triplet_form_append_element(t,k,p,q);
		else
			triplet_form_append_element(t,k,q,p);
		(*size)++;
		cur = cur->next;
	}
  
	if (*size == nedges)
		error = 0;
	goto final;

error_malloc:
	fprintf(stderr,"Error in a malloc or calloc call\n");
	exit(EXIT_FAILURE);

 final:
  if (error == 1)
    free_vector(edges);
  return edges;
}

static PTriangle parseTriangles(xmlNodePtr node,PNode nodes,int nnodes,PEdge edges,SparseMatrix sm,int *size)
{
	PTriangle triangles;
	int found, error;
	size_t ntriangles;
	xmlNodePtr cur;

	error = 1;
	triangles = NULL;
	ntriangles = get_int_value(node,"size",&found,XML_VALUE_REQUIRED);
	make_vector(triangles,ntriangles);
	*size = 0;

	cur = node->xmlChildrenNode;
	while (cur != NULL)
	{
		int index, p, q, r, i;

		if (! AREEQUAL(cur->name,"triangle"))
		{
			cur = cur->next;
			continue;
		}
		index = get_int_value(cur,"index",&found,XML_VALUE_REQUIRED);
		p = get_int_value(cur,"v0",&found,XML_VALUE_REQUIRED);
		q = get_int_value(cur,"v1",&found,XML_VALUE_REQUIRED);
		r = get_int_value(cur,"v2",&found,XML_VALUE_REQUIRED);

		if ((index >= ntriangles) || (p >= nnodes) || (q >= nnodes) || (r >= nnodes))
			goto final;
		if ((p == q) || (p == r) || (q == r))
			goto final;

		PTriangle tr = &triangles[index];
		tr->triangleno = index;
		tr->n[0] = &nodes[p];
		tr->n[1] = &nodes[q];
		tr->n[2] = &nodes[r];
		for (i = 0;i < 3;i++)
		{
			int j = (i + 1) % 3;
			int k = (i + 2) % 3;

			int n1 = tr->n[j]->nodeno;
			int n2 = tr->n[k]->nodeno;
			tr->e[i] = EDGE(n1,n2);
			tr->ex[i] = tr->n[k]->x - tr->n[j]->x;
			tr->ey[i] = tr->n[k]->y - tr->n[j]->y;
		}
		tr->area = 0.5 * fabs(tr->ex[0] * tr->ey[1] - tr->ex[1] * tr->ey[0]);
		(*size)++;
		cur = cur->next;
	}
  
	if (*size == ntriangles)
		error = 0;
	goto final;

error_malloc:
	fprintf(stderr,"Error in a malloc or calloc call\n");
	exit(EXIT_FAILURE);
  
final:
	if (error == 1)
		free_vector(triangles);
	return triangles;
}

DataMesh parseXMLMeshDocument(char *docname)
{
	xmlDocPtr doc;
	xmlNodePtr cur;
	PNode nodes;
	PEdge edges;
	PTriangle triangles;
	DataMesh mesh;
	int error, nnodes, nedges, ntriangles;
	TripletForm tf;
	SparseMatrix sm;
  
	nodes = NULL;
	edges = NULL;
	triangles = NULL;
	doc = NULL;
	cur = NULL;
	sm = NULL;
	tf = NULL;
	error = 1;
  
	if ((mesh = mesh_new()) == NULL)
		goto final;
  
	if ((doc = xmlParseFile(docname)) == NULL)
		goto final;

	if ((cur = xmlDocGetRootElement(doc)) == NULL)
		goto final;

	if (! AREEQUAL(cur->name,"mesh"))
		goto final;

	cur = cur->xmlChildrenNode;
	while (cur != NULL)
	{
		if (AREEQUAL(cur->name,"vertices"))
		{
			if ((nodes = parseVertices(cur,&nnodes)) == NULL)
				goto final;
		}
		else if (AREEQUAL(cur->name,"edges"))
		{
			if ((tf = triplet_form_init(1024)) == NULL)
				goto final;
			if ((edges = parseEdges(cur,nodes,nnodes,tf,&nedges)) == NULL)
	    		goto final;
		}
		else if (AREEQUAL(cur->name,"triangles"))
		{
			if (tf == NULL)
				goto final;
			if ((sm = sparsem_from_triplet_form(tf,ADDITION)) == NULL)
				goto final;
			if ((triangles = parseTriangles(cur,nodes,nnodes,edges,sm,&ntriangles)) == NULL)
	    		goto final;
	    }
		cur = cur->next;
	}

	mesh->nodes = nodes;
	mesh->edges = edges;
	mesh->triangles = triangles;
	mesh->nnodes = nnodes;
	mesh->nedges = nedges;
	mesh->ntriangles = ntriangles;
	error = 0;

final:
	if (error == 1)
	{
		freeDataMesh(mesh);
		free_vector(nodes);
		free_vector(edges);
		free_vector(triangles);
	}
	freeTripletForm(tf);
	freeSparseMatrix(sm);
	if (doc != NULL)
		xmlFreeDoc(doc);
	return mesh;
}

