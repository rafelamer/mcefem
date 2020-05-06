/**************************************************************************************
* Filename:   mesh.c
* Author:     Rouben Rostamian
*             https://userpages.umbc.edu/~rostamia/cbook/
*
* Modifications by:
*
* Disclaimer: This code is presented "as is" and it has been written to 
*             implement the Finite Element Method in dimension 2. It
*             has been writen for educational purposes.
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
*	          See https://www.gnu.org/licenses/
***************************************************************************************/
#include <tfgfem.h>
#include <triangle.h>
#include <math.h>
#include <time.h>

DataMesh mesh_new()
{
	DataMesh m;
	if ((m = (DataMesh)malloc(sizeof(data_mesh))) == NULL)
		return NULL;
	return m;
}

void free_data_mesh(DataMesh *m)
{
	if (*m == NULL)
		return;
	free_vector((*m)->nodes);
	free_vector((*m)->edges);
	free_vector((*m)->triangles);
	free(*m);
	*m = NULL;
}

static struct triangulateio *region_to_triangle(DataRegion r)
{
	int i;
	struct triangulateio *in = malloc(sizeof *in);
	if (in == NULL)
		goto error_malloc;

	/* 
		Process points 
	*/
	in->numberofpoints = r->npoints;
	make_vector(in->pointlist, 2 * in->numberofpoints);
	for (i = 0; i < in->numberofpoints; i++)
	{
		in->pointlist[2*i]   = r->points[i].x;
		in->pointlist[2*i+1] = r->points[i].y;
	}

	make_vector(in->pointmarkerlist, in->numberofpoints);
	for (i = 0; i < in->numberofpoints; i++)
		in->pointmarkerlist[i] = r->points[i].bc;

	in->numberofpointattributes = 0;
  
	/* 
		Process segments 
	*/
	in->numberofsegments = r->nsegments;
	make_vector(in->segmentlist, 2 * in->numberofsegments);
	for (i = 0; i < in->numberofsegments; i++)
	{
		in->segmentlist[2*i]   = r->segments[i].point_no_1;
		in->segmentlist[2*i+1] = r->segments[i].point_no_2;
	}

	make_vector(in->segmentmarkerlist, in->numberofsegments);
	for (i = 0; i < in->numberofsegments; i++)
		in->segmentmarkerlist[i] = r->segments[i].bc;

	/*
		Process holes 
	*/
  	in->numberofholes = r->nholes;
	if (in->numberofholes != 0)
	{
		make_vector(in->holelist, 2 * in->numberofholes);
		for (i = 0; i < in->numberofholes; i++)
		{
			in->holelist[2*i]   = r->holes[i].x;
			in->holelist[2*i+1] = r->holes[i].y;
		}
    }
	else
	{
		in->holelist = NULL;
	}
	in->numberofregions = 0;

	return in;

error_malloc:
	fprintf(stderr,"Error in a malloc or calloc call\n");
	exit(EXIT_FAILURE);
}

static struct triangulateio *do_triangulate(struct triangulateio *in, double a)
{
	char *opts_string = "Qzpeq30a%f";
	char opts[64];
	struct triangulateio *out = malloc(sizeof *out);
	if (out == NULL)
		goto error_malloc;

	out->pointlist = NULL;
	out->pointmarkerlist = NULL;
	out->edgelist = NULL;
	out->edgemarkerlist = NULL;
	out->trianglelist = NULL;
 	out->segmentlist = NULL;
 	out->segmentmarkerlist = NULL;
	sprintf(opts, opts_string, a);
	triangulate(opts, in, out, NULL);
 	return out;

error_malloc:
	fprintf(stderr,"Error in a malloc or calloc call\n");
	exit(EXIT_FAILURE);
}

/* 
	This is ugly. 
	There ought to be a better way but I don't see it 
*/
static void assign_elem_edges(PTriangle triangles,int ntriangles,edge *edges,int nedges,DOUBLE *area)
{
	int i, r, s;
	*area = 0.0;
	for (r = 0; r < ntriangles; r++)
	{
		PTriangle tr = &triangles[r];
		for (i = 0; i < 3; i++)
		{
			int j = i % 3;
			int k = (i + 1) % 3;
			int n1 = tr->n[j]->nodeno;
			int n2 = tr->n[k]->nodeno;
			tr->ex[i] = tr->n[k]->x - tr->n[j]->x;
			tr->ey[i] = tr->n[k]->y - tr->n[j]->y;
			for (s = 0; s < nedges; s++)
			{
				int m1 = edges[s].n[0]->nodeno;
				int m2 = edges[s].n[1]->nodeno;
				if ((m1 == n1 && m2 == n2) || (m1 == n2 && m2 == n1))
				{
					tr->e[i] = &edges[s];
					break;
				}
			}
		}
		tr->area = 0.5 * fabs(tr->ex[0] * tr->ey[1] - tr->ex[1] * tr->ey[0]);
		*area += tr->area;
	}
}

static void test_bc_type(int bc,DataMesh mesh)
{
	switch (bc)
	{
		case FEM_BC_DIRICHLET:
			mesh->hasDirichletBC = 1;
			break;
		case FEM_BC_NEUMANN:
			mesh->hasNeumannBC = 1;
			break;
		case FEM_BC_ROBIN:
			mesh->hasRobinBC = 1;
	}
}

static DataMesh triangle_to_mesh(struct triangulateio *out)
{
	node *nodes;
	edge *edges;
	PTriangle triangles;
	int i, nnodes, nedges, ntriangles;
	int nbnodes, nbedges;
	DataMesh mesh;
	if ((mesh =  mesh_new()) == NULL)
		return NULL;

	mesh->hasDirichletBC = 0;
	mesh->hasNeumannBC = 0;
	mesh->hasRobinBC = 0;
	nbnodes = nbedges = 0;
	nnodes = out->numberofpoints;
	make_vector(nodes, nnodes);
	for (i = 0; i < out->numberofpoints; i++)
	{
		nodes[i].nodeno = i;
		nodes[i].x = out->pointlist[2*i];
		nodes[i].y = out->pointlist[2*i+1];
		nodes[i].bc = out->pointmarkerlist[i];
		nodes[i].bnodeno = -1;
		if (nodes[i].bc >= FEM_BC_DIRICHLET)
		{
			nodes[i].bnodeno = nbnodes++;
			test_bc_type(nodes[i].bc,mesh);
		}
	}

	nedges = out->numberofedges;
	make_vector(edges, nedges);
	for (i = 0; i < nedges; i++)
	{
		edges[i].edgeno = i;
		edges[i].n[0] = &nodes[out->edgelist[2*i]];
		edges[i].n[1] = &nodes[out->edgelist[2*i+1]];
		edges[i].bc = out->edgemarkerlist[i];
		edges[i].bedgeno = -1;
		if (edges[i].bc >= FEM_BC_DIRICHLET)
		{
			edges[i].bedgeno = nbedges++;
			test_bc_type(edges[i].bc,mesh);
		}
	}
	ntriangles = out->numberoftriangles;
	make_vector(triangles, ntriangles);
	for (i = 0; i < ntriangles; i++)
	{
		triangles[i].triangleno = i;
		triangles[i].n[0] = &nodes[out->trianglelist[3*i]];
		triangles[i].n[1] = &nodes[out->trianglelist[3*i+1]];
		triangles[i].n[2] = &nodes[out->trianglelist[3*i+2]];
	}
	assign_elem_edges(triangles, ntriangles, edges, nedges, &(mesh->area));

	mesh->nnodes = nnodes;
	mesh->nedges = nedges;
	mesh->ntriangles = ntriangles;
	mesh->nodes = nodes;
	mesh->edges = edges;
	mesh->triangles = triangles;
	mesh->nbnodes = nbnodes;
	mesh->nbedges = nbedges;

	return mesh;

error_malloc:
	fprintf(stderr,"Error in a malloc or calloc call\n");
	exit(EXIT_FAILURE);	
}

static void free_triangle_in_structure(struct triangulateio *in)
{
	free_vector(in->pointlist);
	free_vector(in->pointmarkerlist);
	free_vector(in->segmentlist);
	free_vector(in->segmentmarkerlist);
	free_vector(in->holelist);
	free(in);
}

static void free_triangle_out_structure(struct triangulateio *out)
{
	free(out->pointlist);
	free(out->pointmarkerlist);
	free(out->edgelist);
	free(out->edgemarkerlist);
	free(out->trianglelist);
	free(out->segmentlist);
	free(out->segmentmarkerlist);
	free(out);
}

DataMesh make_mesh(DataRegion r,DOUBLE a)
{
	struct triangulateio *in, *out;
	DataMesh mesh;
	in = region_to_triangle(r);
	out = do_triangulate(in, a);
	mesh = triangle_to_mesh(out);
	free_triangle_in_structure(in);
	free_triangle_out_structure(out);
	return mesh;
}

void mesh_to_asy(DataMesh mesh,int degree,const char *outfile,double unitsize,int labels)
{
	FILE *fp;
	time_t now;
	int i;
	DOUBLE dot;
  
	if ((fp = fopen(outfile, "w")) == NULL)
	{
		fprintf(stderr, "cannot open file %s for writing\n", outfile);
		return;
	}
	fprintf(fp, "/*\n");
	fprintf(fp, " * Title: %s\n", outfile);
	fprintf(fp, " * Creator: Study and implementation in C of the Finite Element Method in dimension 2, %s\n", __FILE__);
	now = time(NULL);
	fprintf(fp, " * CreationDate: %s", ctime(&now));
	fprintf(fp,"*/\n\n");
	fprintf(fp,"usepackage(\"times\");\n");
	fprintf(fp,"usepackage(\"mtpro2\");\n\n");

	fprintf(fp,"unitsize(%.8gcm);\n",unitsize);
	dot = 0.06667 * mesh->area;
	if (dot < 0.2)
		dot = 0.2;
	dot *= unitsize;
	fprintf(fp,"pen y = yellow+linewidth(%.3gpt);\n",dot);
	fprintf(fp,"pen r = red+linewidth(%.3gpt);\n",dot);
	fprintf(fp,"pen b = orange+linewidth(%.3gpt);\n",dot);
	fprintf(fp,"pen g = green+linewidth(%.3gpt);\n",dot);
	
	for (i = 0; i < mesh->ntriangles; i++)
	{
		int j, k, l;
		PTriangle tr = &(mesh->triangles[i]);

		fprintf(fp,"filldraw((%.16g,%.16g)--(%.16g,%.16g)--(%.16g,%.16g)--cycle,fillpen=lightgrey,drawpen=red+linewidth(%.3gpt));\n",
				tr->n[0]->x	,tr->n[0]->y,
				tr->n[1]->x	,tr->n[1]->y,
				tr->n[2]->x	,tr->n[2]->y,
				0.01 * mesh->area);
		
		if (degree < 3)
			continue;

		for(j = 1;j < degree - 1;j++)
		{
			for (k = 1;k < degree - j;k++)
			{
				DOUBLE u[2], v[2];
				u[0] = (tr->n[1]->x - tr->n[0]->x) / degree;
				u[1] = (tr->n[1]->y - tr->n[0]->y) / degree;
				v[0] = (tr->n[2]->x - tr->n[0]->x) / degree;
				v[1] = (tr->n[2]->y - tr->n[0]->y) / degree;
				fprintf(fp,"dot((%.16g,%.16g),y);\n",
						tr->n[0]->x + j * u[0] + k * v[0],
						tr->n[0]->y + j * u[1] + k * v[1]); 
				if (labels != 0)
				{
					int n = global_node_number(mesh,tr,degree,j,k);
					fprintf(fp,"label(\"$u_{%d}$\",(%.16g,%.16g),S);\n",n,
						    tr->n[0]->x + j * u[0] + k * v[0],
						    tr->n[0]->y + j * u[1] + k * v[1]);
				}
			}
		}
	}
	
	for(i = 0;i < mesh->nedges;i++)
	{
		PEdge eg = &(mesh->edges[i]);
		if(eg->bc < FEM_BC_DIRICHLET)
			continue;
		fprintf(fp,"draw((%.16g,%.16g)--(%.16g,%.16g),blue+linewidth(%.3gpt));\n",
			eg->n[0]->x,
			eg->n[0]->y,
			eg->n[1]->x,
			eg->n[1]->y,
			0.025 * mesh->area);
	}
	
	if (degree < 1)
		goto final;
	
	for(i = 0;i < mesh->nnodes;i++)
	{
		PNode nd = &(mesh->nodes[i]);
		if (nd->bc == FEM_BC_DIRICHLET)
			fprintf(fp,"dot((%.16g,%.16g),b);\n",nd->x,nd->y);
		else if (nd->bc == FEM_BC_NEUMANN)
			fprintf(fp,"dot((%.16g,%.16g),r);\n",nd->x,nd->y);
		else if (nd->bc == FEM_BC_ROBIN)
			fprintf(fp,"dot((%.16g,%.16g),g);\n",nd->x,nd->y);		
		else
			fprintf(fp,"dot((%.16g,%.16g),y);\n",nd->x,nd->y);
		if (labels != 0)
		{
			if (nd->bc < FEM_BC_DIRICHLET)
				fprintf(fp,"label(\"$u_{%d}$\",(%.16g,%.16g),1.2*S);\n",nd->nodeno,nd->x,nd->y);
			else
				fprintf(fp,"label(\"$u_{%d},\\dot u_{%d}$\",(%.16g,%.16g),1.2*S);\n",nd->nodeno,nd->bnodeno,nd->x,nd->y);
		}
	}
	
	if (degree < 2)
		goto final;
	
	for(i = 0;i < mesh->nedges;i++)
	{
		DOUBLE v0, v1;
		int k;
		PEdge eg = &(mesh->edges[i]);
		
		v0 = (eg->n[1]->x - eg->n[0]->x) / degree;
		v1 = (eg->n[1]->y - eg->n[0]->y) / degree;
		for (k = 1; k < degree;k++)
		{
			if (eg->bc == FEM_BC_DIRICHLET)
				fprintf(fp,"dot((%.16g,%.16g),b);\n",eg->n[0]->x + k * v0,eg->n[0]->y + k * v1);
			else if (eg->bc == FEM_BC_NEUMANN)
				fprintf(fp,"dot((%.16g,%.16g),r);\n",eg->n[0]->x + k * v0,eg->n[0]->y + k * v1);
			else if (eg->bc == FEM_BC_ROBIN)
				fprintf(fp,"dot((%.16g,%.16g),g);\n",eg->n[0]->x + k * v0,eg->n[0]->y + k * v1);
			else
				fprintf(fp,"dot((%.16g,%.16g),y);\n",eg->n[0]->x + k * v0,eg->n[0]->y + k * v1);
			if (labels != 0)
			{	int n = global_node_number_at_edge(mesh,eg,degree,k);
				if (eg->bc < FEM_BC_DIRICHLET)
				{
					int n = global_node_number_at_edge(mesh,eg,degree,k);
					fprintf(fp,"label(\"$u_{%d}$\",(%.16g,%.16g),1.2*S);\n",n,eg->n[0]->x + k * v0,eg->n[0]->y + k * v1);
				}
				else
				{
					int m = boundary_node_number_at_edge(mesh,eg,degree,k);
					fprintf(fp,"label(\"$u_{%d},\\dot u_{%d}$\",(%.16g,%.16g),1.2*S);\n",n,m,eg->n[0]->x + k * v0,eg->n[0]->y + k * v1);	
				}
			}
		}
	}
	
final:
  fclose(fp);
}

void mesh_to_asy_boundary(DataMesh mesh,const char *outfile)
{
	FILE *fp;
	time_t now;
	int i;
  
	if ((fp = fopen(outfile, "w")) == NULL)
	{
		fprintf(stderr, "cannot open file %s for writing\n", outfile);
		return;
	}
	fprintf(fp, "/*\n");
	fprintf(fp, " * Title: (%s)\n", outfile);
	fprintf(fp, " * Creator: Study and implementation in C of the Finite Element Method in dimension 2, %s\n", __FILE__);
	now = time(NULL);
	fprintf(fp, " * CreationDate: %s", ctime(&now));
	fprintf(fp,"*/\n");
	fprintf(fp,"unitsize(1.0cm);\n\n");
	fprintf(fp,"arrowbar a = Arrow(size=1mm);\n");

	for (i = 0; i < mesh->ntriangles; i++)
	{
		int j;
		PTriangle tr = &(mesh->triangles[i]);

		fprintf(fp,"filldraw((%.16g,%.16g)--(%.16g,%.16g)--(%.16g,%.16g)--cycle,fillpen=lightgrey,drawpen=red+linewidth(0.1pt));\n",
				tr->n[0]->x	,tr->n[0]->y,
				tr->n[1]->x	,tr->n[1]->y,
				tr->n[2]->x	,tr->n[2]->y);

		for (j = 0;j < 3;j++)
		{
			PEdge eg;
			eg = tr->e[j];
			if (eg->bc >= 2)
				fprintf(fp,"draw((%.16g,%.16g)--(%.16g,%.16g),blue+linewidth(0.1666pt),a);\n",
						eg->n[0]->x,eg->n[0]->y,eg->n[1]->x,eg->n[1]->y);
		}

	}
	
final:
	fclose(fp);
}

void mesh_to_asy_region(DataMesh mesh,const char *outfile)
{
	FILE *fp;
	time_t now;
	int i;
  
	if ((fp = fopen(outfile, "w")) == NULL)
	{
		fprintf(stderr, "cannot open file %s for writing\n", outfile);
		return;
	}
	fprintf(fp, "/*\n");
	fprintf(fp, " * Title: (%s)\n", outfile);
	fprintf(fp, " * Creator: Study and implementation in C of the Finite Element Method in dimension 2, %s\n", __FILE__);
	now = time(NULL);
	fprintf(fp, " * CreationDate: %s", ctime(&now));
	fprintf(fp,"*/\n");
	fprintf(fp,"unitsize(1.0cm);\n\n");
	fprintf(fp,"pen g = lightgrey;\n");
	fprintf(fp,"pen f = lightgrey+linewidth(0pt);\n");
	fprintf(fp,"pen b = blue+linewidth(0.2666pt);\n\n");

	for (i = 0; i < mesh->ntriangles; i++)
	{
		int j;
		PTriangle tr = &(mesh->triangles[i]);

		fprintf(fp,"filldraw((%.16g,%.16g)--(%.16g,%.16g)--(%.16g,%.16g)--cycle,fillpen=g,drawpen=f);\n",
				tr->n[0]->x	,tr->n[0]->y,
				tr->n[1]->x	,tr->n[1]->y,
				tr->n[2]->x	,tr->n[2]->y);
	}
	for (i = 0; i < mesh->ntriangles; i++)
	{
		int j;
		PTriangle tr = &(mesh->triangles[i]);
		for (j = 0;j < 3;j++)
		{
			PEdge eg;
			eg = tr->e[j];
			if (eg->bc >= 2)
				fprintf(fp,"draw((%.16g,%.16g)--(%.16g,%.16g),b);\n",
						eg->n[0]->x,eg->n[0]->y,eg->n[1]->x,eg->n[1]->y);
		}
	}
	
final:
	fclose(fp);
}

void print_triangle(DataMesh mesh,int n) 
{
	PTriangle tr = &(mesh->triangles[n]);
	printf("Triangle number %d\n",n);
	printf("   Node 0: Number: %d   (%.16g,%.16g)\n",tr->n[0]->nodeno,tr->n[0]->x,tr->n[0]->y);
	printf("   Node 1: Number: %d   (%.16g,%.16g)\n",tr->n[1]->nodeno,tr->n[1]->x,tr->n[1]->y);
	printf("   Node 2: Number: %d   (%.16g,%.16g)\n\n",tr->n[2]->nodeno,tr->n[2]->x,tr->n[2]->y);
	printf("   Edge 0: Number: %d Begin: %d  End: %d\n",tr->e[0]->edgeno,tr->e[0]->n[0]->nodeno,tr->e[0]->n[1]->nodeno);
	printf("   Edge 1: Number: %d Begin: %d  End: %d\n",tr->e[1]->edgeno,tr->e[1]->n[0]->nodeno,tr->e[1]->n[1]->nodeno);
	printf("   Edge 2: Number: %d Begin: %d  End: %d\n\n",tr->e[2]->edgeno,tr->e[2]->n[0]->nodeno,tr->e[2]->n[1]->nodeno);
}

void print_list_of_distinct_points(DataMesh mesh,int degree,const char *outfile)
{
	int i, count;
	FILE *fp;

	if (degree < 1)
	{
		fprintf(stderr, "Degree must be a positive integer\n");
		return;
	}
	
	if ((fp = fopen(outfile, "w")) == NULL)
	{
		fprintf(stderr, "cannot open file %s for writing\n", outfile);
		return;
	}

	count = 0;
	/*
		Nodes
	*/
	for(i = 0;i < mesh->nnodes;i++)
	{
		PNode nd = &(mesh->nodes[i]);
		fprintf(fp,"%d  Node: %d (%.16g,%.16g)\n",count++,nd->nodeno,nd->x,nd->y);
		if (nd->bnodeno >= 0)
			fprintf(fp,"          Boundary node: %d\n",nd->bnodeno);
	}

	if (degree < 2)
		goto final;
	/*
		Interior points of the edges
	*/
	for(i = 0;i < mesh->nedges;i++)
	{
		DOUBLE v0, v1;
		int k;
		PEdge eg = &(mesh->edges[i]);
		
		v0 = (eg->n[1]->x - eg->n[0]->x) / degree;
		v1 = (eg->n[1]->y - eg->n[0]->y) / degree;
		for (k = 1; k < degree;k++)
		{ 
			fprintf(fp,"%d  Edge: %d k=%d (%.16g,%.16g)\n",count++,eg->edgeno,k,eg->n[0]->x + k * v0,eg->n[0]->y + k * v1);
			if (eg->bedgeno >= 0)
				fprintf(fp,"          Boundary node: %d\n",mesh->nbnodes + (degree - 1) * eg->bedgeno + k - 1);
		}	
	}

	if (degree < 3)
		goto final;
	/*
		Interior points of the triangle
	*/
	for (i = 0; i < mesh->ntriangles; i++)
	{
		int j, k;
		PTriangle tr = &(mesh->triangles[i]);
		for(j = 1;j < degree - 1;j++)
		{
			for (k = 1;k < degree - j;k++)
			{
				DOUBLE u[2], v[2];
				u[0] = (tr->n[1]->x - tr->n[0]->x) / degree;
				u[1] = (tr->n[1]->y - tr->n[0]->y) / degree;
				v[0] = (tr->n[2]->x - tr->n[0]->x) / degree;
				v[1] = (tr->n[2]->y - tr->n[0]->y) / degree;
				fprintf(fp,"%d  Triangle: %d  Local coordinates: (%d,%d)  Global coordinates: (%.16g,%.16g)\n",
						count++,tr->triangleno,j,k,
						tr->n[0]->x + j * u[0] + k * v[0],
						tr->n[0]->y + j * u[1] + k * v[1]);
			}
		}
	}
	
final:
  fclose(fp);
}

int global_node_number(DataMesh mesh,PTriangle tr,int degree,int m,int n) 
{
	int nodes, edges;

	if ((m < 0) || (n < 0))
		return -1;
	if ((m > degree) || (n > degree - m))
		return -1;
	/*
		It's a primary node point
	*/
	if ((m == 0) && (n == 0))
		return tr->n[0]->nodeno;
	if ((m == degree) && (n == 0))
		return tr->n[1]->nodeno;
	if ((m == 0) && (n == degree))
		return tr->n[2]->nodeno;
	/*
		It's an interior edge point
	*/
	nodes = mesh->nnodes;
	PEdge eg = NULL;
	if (m == 0)
	{
		eg = tr->e[2];
		if (tr->n[0]->nodeno == eg->n[0]->nodeno)
			return nodes + (degree - 1) * eg->edgeno + n - 1;
		else
			return nodes + (degree - 1) * eg->edgeno + (degree - n) - 1;
	}
	else if (n == 0)
	{
		eg = tr->e[0];
		if (tr->n[0]->nodeno == eg->n[0]->nodeno)
			return nodes + (degree - 1) * eg->edgeno + m - 1;
		else
			return nodes + (degree - 1) * eg->edgeno + (degree - m) - 1;
	}
	else if (m + n == degree)
	{
		eg = tr->e[1];
		if (tr->n[1]->nodeno == eg->n[0]->nodeno)
			return nodes + (degree - 1) * eg->edgeno + n - 1;
		else
			return nodes + (degree - 1) * eg->edgeno + (degree - n) - 1;
	}
	/*
		It's an interior points of the triangle
	*/
	nodes = mesh->nnodes + (degree - 1) * mesh->nedges + tr->triangleno * (degree - 2) * (degree - 1) / 2; 
	return nodes + ((m - 1) * (2 * (degree - 1) - m))/2 + n - 1;
}

int global_node_number_at_edge(DataMesh mesh,PEdge edge,int degree,int m)
{
	if ((m < 0) || (m > degree))
		return -1;
	/*
		It's a primary node point
	*/
	if (m == 0)
		return edge->n[0]->nodeno;
	if (m == degree)
		return edge->n[1]->nodeno;

	/*
		It's an interior boundary edge point
	*/
	return mesh->nnodes + (degree - 1) * edge->edgeno + m - 1;
} 

int boundary_node_number(DataMesh mesh,PTriangle tr,int degree,int m,int n)
{
	int bnodes;
	if ((m < 0) || (n < 0))
		return -1;
	if ((m > degree) || (n > degree - m))
		return -1;
	
	/*
		It's a primary node point
	*/
	if ((m == 0) && (n == 0))
		return tr->n[0]->bnodeno;
	if ((m == degree) && (n == 0))
		return tr->n[1]->bnodeno;
	if ((m == 0) && (n == degree))
		return tr->n[2]->bnodeno;

	/*
		It's an interior boundary edge point
	*/
	bnodes = mesh->nbnodes;
	if (n == 0)
	{
		PEdge eg;
		eg = tr->e[0];	
		if (eg->bedgeno >= 0)
			return bnodes + (degree - 1) * eg->bedgeno + m - 1;
	}
	if (m == 0)
	{
		PEdge eg;
		eg = tr->e[1];
		if (eg->bedgeno >= 0)	
			return bnodes + (degree - 1) * eg->bedgeno + n - 1;
	}
	if (m + n == degree)
	{
		PEdge eg;
		eg = tr->e[2];	
		if (eg->bedgeno >= 0)
			return bnodes + (degree - 1) * eg->edgeno + m - 1;
	}
	/*
		It's not a bounadry node
	*/
	return -1;
}

int boundary_node_number_at_edge(DataMesh mesh,PEdge edge,int degree,int m)
{	
	if ((m < 0) || (m > degree))
		return -1;
	/*
		It's a primary node point
	*/
	if (m == 0)
		return edge->n[0]->bnodeno;
	if (m == degree)
		return edge->n[1]->bnodeno;
	/*
		It's an interior boundary edge point
	*/
	return mesh->nbnodes + (degree - 1) * edge->bedgeno + m - 1;
}

int boundary_condition_at_node_edge(PEdge edge,int degree,int m)
{	
	if ((m < 0) || (m > degree))
		return -1;
	/*
		It's a primary node point
	*/
	if (m == 0)
		return edge->n[0]->bc;
	if (m == degree)
		return edge->n[1]->bc;
	/*
		It's an interior boundary edge point
	*/
	return edge->bc;
}

void node_coordinates(PTriangle tr,int degree,int m,int n,PNode r) 
{
	PNode nd = NULL;

	if ((m < 0) || (n < 0))
	{
		r->x = NAN;
		r->y = NAN;
		return;
	}
	if ((m > degree) || (n > degree - m))
	{
		r->x = NAN;
		r->y = NAN;
		return;
	}
	PEdge eg = NULL;
	DOUBLE u[2], v[2];
	u[0] = (tr->n[1]->x - tr->n[0]->x) / degree;
	u[1] = (tr->n[1]->y - tr->n[0]->y) / degree;
	v[0] = (tr->n[2]->x - tr->n[0]->x) / degree;
	v[1] = (tr->n[2]->y - tr->n[0]->y) / degree;
	r->x = tr->n[0]->x + m * u[0] + n * v[0];
	r->y = tr->n[0]->y + m * u[1] + n * v[1];
}

void node_coordinates_at_edge(PEdge edge,int degree,int m,PNode r)
{
	if ((m < 0) || (m > degree))
	{
		r->x = NAN;
		r->y = NAN;
		return;
	}
	r->x = edge->n[0]->x + m * (edge->n[1]->x - edge->n[0]->x) / degree;
	r->y = edge->n[0]->y + m * (edge->n[1]->y - edge->n[0]->y) / degree;
}


void print_list_of_points_per_triangle(DataMesh mesh,int degree,const char *outfile)
{
	int tn;
	FILE *fp;

	if (degree < 1) 
	{
		fprintf(stderr, "Degree must be a positive integer\n");
		return;
	}
	
	if ((fp = fopen(outfile, "w")) == NULL)
	{
		fprintf(stderr, "cannot open file %s for writing\n", outfile);
		return;
	}
	for (tn = 0; tn < mesh->ntriangles; tn++)
	{
		int i, j;
		PTriangle tr = &(mesh->triangles[tn]);
		fprintf(fp,"Triangle number %d\n",tn);
		node nd;

		for(i = 0;i <= degree;i++)
		{
			for (j = 0;j <= degree - i;j++)
			{

				int nodeno = global_node_number(mesh,tr,degree,i,j);
				node_coordinates(tr,degree,i,j,&nd);
				fprintf(fp,"   Local coordinates: (%d,%d)   Unknown: %d  Global coordinates: (%.16g,%.16g)\n",i,j,nodeno,nd.x,nd.y);
			}
		}
	}
}

void print_mesh_data(DataMesh mesh,int degree) {
	int i, j;
	printf("Number of points: %d\n",mesh->nnodes);
	printf("Number of edges: %d\n",mesh->nedges);
	printf("Number of triangles: %d\n\n",mesh->ntriangles);
	printf("Number of boundary nodes: %d\n\n",mesh->nbnodes);
	printf("Number of boundary edges: %d\n\n",mesh->nbedges);
	printf("Number of primary unknowns: %d\n",number_of_u_unknowns(mesh,degree));
	printf("Number of seconary unknowns: %d\n",number_of_dotu_unknowns(mesh,degree));
	printf("Boundary nodes\n");
	for (i = 0;i < mesh->nnodes;i++)
	{
		PNode nd;
		nd = &(mesh->nodes[i]);
		if (nd->bnodeno >= 0)
			printf("    Boundary node number: %d\n",nd->bnodeno);
	}
	printf("Boundary nodes at edges\n");
	for (i = 0;i < mesh->nedges;i++)
	{
		PEdge eg;
		eg = &(mesh->edges[i]);
		if (eg->bedgeno >= 0) 
		{
			printf("    Boundary edge number: %d\n",eg->bedgeno);
			for (j = 1;j < degree;j++)
				printf("        Boundary node in the edge number: %d\n",
						boundary_node_number_at_edge(mesh,eg,degree,j));
		}
	}
}

int number_of_u_unknowns(DataMesh mesh,int degree)
{
	int primary;

	primary = mesh->nnodes;

	if (degree <= 1)
		goto final;

	if (degree == 2)
	{
		primary += mesh->nedges;
		goto final;
	}
	primary += (degree - 1) * mesh->nedges;
	primary += (degree - 2) * (degree - 1) * mesh->ntriangles / 2;

final:
	return primary;
}

int number_of_dotu_unknowns(DataMesh mesh,int degree)
{
	int secondary;

	secondary = mesh->nbnodes;

	if (degree <= 1)
		goto final;

	if (degree == 2)
	{
		secondary += mesh->nbedges;
		goto final;
	}
	secondary += (degree - 1) * mesh->nbedges;

final:
	return secondary;
}

int number_of_unknowns(DataMesh mesh,int degree)
{
	return number_of_u_unknowns(mesh,degree) + number_of_dotu_unknowns(mesh,degree); 
}
