/**************************************************************************************
* Filename:   fem2.c
* Authors:     
* Copyright:  
* Disclaimer: This code is presented "as is" and it has been written to 
*             implement the Finite Element Method in dimensions 1 and 2.
*             It has been writen educational purposes.
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
#include <math.h>

#define Value(m,n) (lv2[(m)][(n)][count].value)
#define DXValue(m,n) (a22 * lv2[(m)][(n)][count].dxvalue - a21 * lv2[(m)][(n)][count].dyvalue)
#define DYValue(m,n) (-a12 * lv2[(m)][(n)][count].dxvalue + a11 * lv2[(m)][(n)][count].dyvalue)

#define load_vector_over_one_triangle(m1,n1)  {                                        \
	DOUBLE x[3], y[3];                                                                 \
	DOUBLE la1, la2, X, Y, a11, a12, a21, a22, det, factor;                            \
	int z, count;                                                                      \
	TWBGaussQuadrature q2;                                                             \
	r = 0.0;                                                                           \
	q2 = spec->qdat2d;                                                                 \
	for (z = 0; z < 3; z++)                                                            \
	{                                                                                  \
		x[z] = tr->n[z]->x;                                                            \
		y[z] = tr->n[z]->y;                                                            \
	}                                                                                  \
	a11 = x[1] - x[0];                                                                 \
	a12 = x[2] - x[0];                                                                 \
	a21 = y[1] - y[0];                                                                 \
	a22 = y[2] - y[0];                                                                 \
	det = a11 * a22 - a12 * a21;                                                       \
	count = 0;                                                                         \
	while (q2->weight != -1)                                                           \
	{                                                                                  \
		DOUBLE t;                                                                      \
		la1 = q2->lambda1;                                                             \
		la2 = q2->lambda2;                                                             \
		X = x[0] + (x[1] - x[0]) * la1 + (x[2] - x[0]) * la2;                          \
		Y = y[0] + (y[1] - y[0]) * la1 + (y[2] - y[0]) * la2;                          \
		t = q2->weight * runFunctionTable2D(f,X,Y) * Value(m1,n1);                     \
		r += t;                                                                        \
		q2++;                                                                          \
		count++;                                                                       \
	}                                                                                  \
	r *= 0.5 * tr->area;                                                               \
}

#define stiffnes_matrix_over_one_triangle(m1,n1,m2,n2) {                                   \
	DOUBLE x[3], y[3];                                                                     \
	DOUBLE la1, la2, X, Y, sum, a11, a12, a21, a22, det, factor;                           \
	int z, count;                                                                          \
	TWBGaussQuadrature q2;                                                                 \
	r1 = r2 = r3 = r4 = 0.0;                                                               \
	q2 = spec->qdat2d;                                                                     \
	for (z = 0; z < 3; z++)                                                                \
	{                                                                                      \
		x[z] = tr->n[z]->x;                                                                \
		y[z] = tr->n[z]->y;                                                                \
	}                                                                                      \
	a11 = x[1] - x[0];                                                                     \
	a12 = x[2] - x[0];                                                                     \
	a21 = y[1] - y[0];                                                                     \
	a22 = y[2] - y[0];                                                                     \
	det = a11 * a22 - a12 * a21;                                                           \
	count = 0;                                                                             \
	while (q2->weight != -1)                                                               \
	{                                                                                      \
		DOUBLE t;                                                                          \
		la1 = q2->lambda1;                                                                 \
		la2 = q2->lambda2;                                                                 \
		X = x[0] + (x[1] - x[0]) * la1 + (x[2] - x[0]) * la2;                              \
		Y = y[0] + (y[1] - y[0]) * la1 + (y[2] - y[0]) * la2;                              \
		t = q2->weight * runFunctionTable2D(a,X,Y);                                        \
		t *= DXValue(m1,n1) * DXValue(m2,n2) + DYValue(m1,n1) * DYValue(m2,n2);            \
		r1 += t;                                                                           \
		t = q2->weight * runFunctionTable2D(c,X,Y) * Value(m1,n1) * Value(m2,n2);          \
		r3 += t;                                                                           \
		t = q2->weight * runFunctionTable2D(b1,X,Y) * Value(m1,n1) * DXValue(m2,n2);       \
		r2 += t;                                                                           \
		t = q2->weight * runFunctionTable2D(b2,X,Y) * Value(m1,n1) * DYValue(m2,n2);       \
		r2 += t;                                                                           \
        if ((m1 != m2) || (n1 != n2))                                                      \
        {                                                                                  \
			t = q2->weight * runFunctionTable2D(b1,X,Y) * Value(m2,n2) * DXValue(m1,n1);   \
			r4 += t;                                                                       \
			t = q2->weight * runFunctionTable2D(b2,X,Y) * Value(m2,n2) * DYValue(m1,n1);   \
			r4 += t;                                                                       \
		}                                                                                  \
		q2++;                                                                              \
		count++;                                                                           \
	}                                                                                      \
	r1 *= 0.125/tr->area;                                                                  \
    r2 *= 0.25 * ((det > 0) - (det < 0));                                                  \
    r3 *= 0.5 * tr->area;                                                                  \
    r4 *= 0.25 * ((det > 0) - (det < 0));                                                  \
}

#define integral_over_one_edge_two_lagrange(item,edge,m1,n1,r) {     \
	DOUBLE X, Y;                                                     \
	int count = 0;                                                   \
	r = 0.0;                                                         \
	GaussQuadrature q1 = spec->qdat1d;                               \
	DOUBLE v0, v1;                                                   \
	v0 = edge->n[1]->x - edge->n[0]->x;                              \
	v1 = edge->n[1]->y - edge->n[0]->y;                              \
	while(q1->w != -1)                                               \
	{                                                                \
		DOUBLE t;                                                    \
		t = 0.5 * q1->x + 0.5;                                       \
		X = edge->n[0]->x + t * v0;                                  \
		Y = edge->n[0]->y + t * v1;                                  \
		r += q1->w * runFunctionTable2D(item,X,Y) *                  \
					 lv1[m1][count].value * lv1[n1][count].value;    \
		q1++;                                                        \
		count++;                                                     \
	}                                                                \
	r *= 0.5 * sqrt(v0 * v0 + v1 * v1);                              \
}

SystemOfEquations StiffnessMatrixAndLoadVector2D(Specification2D spec)
{
	SystemOfEquations s = NULL;
	TripletForm t = NULL;
	Lagrange1DValues lv1;
	Lagrange2DValues lv2;
	DataMesh mesh;
	int error = 1;
	int us, dotus,tn, degree;
	FunctionTable a, b1, b2, c, f, d, n, A, B, C;

	/*
		Usefull variables
	*/ 
	mesh = spec->mesh;
	lv1 = spec->lv1d;
	lv2 = spec->lv2d;
	degree = spec->degree;
	
	/*
		Table functions
	*/
	HASH_FIND_STR(spec->funs,"a",a);
	HASH_FIND_STR(spec->funs,"b1",b1);
	HASH_FIND_STR(spec->funs,"b2",b2);
	HASH_FIND_STR(spec->funs,"c",c);
	HASH_FIND_STR(spec->funs,"f",f);
	if ((a == NULL) || (b1 == NULL) || (b2 == NULL) || (c == NULL) || (f == NULL)) 
		goto error_malloc;

	HASH_FIND_STR(spec->funs,"d",d);
	HASH_FIND_STR(spec->funs,"n",n);
	HASH_FIND_STR(spec->funs,"A",A);
	HASH_FIND_STR(spec->funs,"B",B);
	HASH_FIND_STR(spec->funs,"C",C);

	if ((mesh->hasDirichletBC != 0) && (d == NULL))
		goto error_malloc;
	if ((mesh->hasNeumannBC != 0) && (n == NULL))
		goto error_malloc;
	if (mesh->hasRobinBC != 0)
		if ((A == NULL) || (B == NULL) || (C == NULL))
			goto error_malloc;	

	/*
		Number of unknowns, system of equations and triplet form
	*/
	us =  number_of_u_unknowns(mesh,degree);
	dotus =  number_of_dotu_unknowns(mesh,degree);
	
	if ((s = (SystemOfEquations)malloc(sizeof(system_of_equations))) == NULL)
		goto error_malloc;

	if ((t = triplet_form_init(SIZEINCREMENT)) == NULL)
		goto error_malloc;

	make_vector(s->F,us + dotus);

	/*
		We start to compute integrals
	*/
	for (tn = 0; tn < mesh->ntriangles; tn++)
	{
		int i, j, k, l;
		/*
			For every pair (i,j), we have a node point. It can be
				1. A node point of the triangulation
				2. An interior edge point
				3. An interior point of the triangle

		*/
		PTriangle tr = &(mesh->triangles[tn]);
		for(i = 0;i <= degree;i++)
		{	
			for (j = 0;j <= degree - i;j++)
			{
				int nij = global_node_number(mesh,tr,degree,i,j);
				DOUBLE r;
				/*
					Load vector: \int_{tr} f(x,y) \phi_{ij}(x,y) dxdy
				*/
				load_vector_over_one_triangle(i,j);
				s->F[nij] += r;
				/*
					For every pair (k,l), we have another node point and we have
					to compute integrals over a triangle with two Lagrange functions
				*/
				for(k = 0;k <= i;k++)
				{
					for (l = 0;l <= degree - k;l++)
					{
						if ((k == i) && (l > j))
							break;
						int nkl = global_node_number(mesh,tr,degree,k,l);
						/*
							r1 = \int_{tr} a(x,y) (\nabla\phi_{ij}(x,y) \cdot \nabla\phi_{kl}(x,y)) dxdy
							r3 = \int_{tr} c(x,y) \phi_{ij}(x,y) \phi_{kl}(x,y) dxdy
							r2 = \int_{tr} (b1(x,y) \frac{\partial}{\partial x}\phi_{kl}(x,y) 
							                + b2(x,y) \frac{\partial}{\partial y}\phi_{kl}(x,y)) \phi_{ij}(x,y) dxdy 
							r4 = \int_{tr} (b1(x,y) \frac{\partial}{\partial x}\phi_{ij}(x,y) 
							                + b2(x,y) \frac{\partial}{\partial y}\phi_{ij}(x,y)) \phi_{kl}(x,y) dxdy 
						*/
						DOUBLE r1, r2, r3, r4;
						stiffnes_matrix_over_one_triangle(i,j,k,l);
						if ((i == k) && (j == l))
						{
							if (r1 + r3 != 0.0)	
								if (! triplet_form_append_element(t,r1+r3,nij,nij))
									goto error_malloc;
							if (r2 != 0)
								if (! triplet_form_append_element(t,r2,nij,nij))
									goto error_malloc;
						}
						else
						{
							if (r1 + r3 != 0.0)
							{	
								if (! triplet_form_append_element(t,r1+r3,nij,nkl))
									goto error_malloc;
								if (! triplet_form_append_element(t,r1+r3,nkl,nij))
									goto error_malloc;
							}
							if (r2 != 0.0)
								if (! triplet_form_append_element(t,r2,nij,nkl))
									goto error_malloc;
							if (r4 != 0.0)
								if (! triplet_form_append_element(t,r4,nkl,nij))
									goto error_malloc;
						}
					}
				}
			}
		}
		/*
			For the edges of the triangle that are boundary edges, we have to compute
			the terms \int_{\Gamma} a(x,y) \phi_{i}(x,y) \partial_{n}u(x,y) ds
			that apperars at equation i
		*/
		for (k = 0;k < 3;k++)
		{
			PEdge edge = tr->e[k];
			if(edge->bedgeno < 0)
				continue;
			for (i = 0;i <= degree;i++)
			{
				int ui = global_node_number_at_edge(mesh,edge,degree,i);
				int dotui = boundary_node_number_at_edge(mesh,edge,degree,i);
				for (j = 0;j <= i;j++)
				{
					DOUBLE r;
					int uj = global_node_number_at_edge(mesh,edge,degree,j);
					int dotuj = boundary_node_number_at_edge(mesh,edge,degree,j);
					integral_over_one_edge_two_lagrange(a,edge,i,j,r);
					if (i == j)
					{
						if (r != 0.0)
							if (! triplet_form_append_element(t,-r,ui,us + dotuj))
								goto error_malloc;
					} 
					else
					{	
						if(r != 0.0) 
						{
							if (! triplet_form_append_element(t,-r,ui,us + dotuj))
								goto error_malloc;
							if (! triplet_form_append_element(t,-r,uj,us + dotui))
								goto error_malloc;
						}
					}
				}
				/*
					Boundary conditions at the points of the edge
				*/
				if (i == degree)
					continue;
				node r;
				int bc = boundary_condition_at_node_edge(edge,degree,i);
				node_coordinates_at_edge(edge,degree,i,&r);
				if (bc == FEM_BC_DIRICHLET)
				{
					if (! triplet_form_append_element(t,1.0,us + dotui,ui))
						goto error_malloc;
					s->F[us + dotui] += runFunctionTable2D(d,r.x,r.y);
				}
				else if (bc == FEM_BC_NEUMANN)
				{
					if (! triplet_form_append_element(t,1.0,us + dotui,us + dotui))
							goto error_malloc;
					s->F[us + dotui] += runFunctionTable2D(n,r.x,r.y);
				}
				else if (bc == FEM_BC_ROBIN)
				{
					if (! triplet_form_append_element(t,runFunctionTable2D(A,r.x,r.y),us + dotui,ui))
							goto error_malloc;
					if (! triplet_form_append_element(t,runFunctionTable2D(B,r.x,r.y),us + dotui,us + dotui))
							goto error_malloc;	
					s->F[us + dotui] += runFunctionTable2D(C,r.x,r.y);
				}
			}
		}	
	}
	/*
		Generate the sparse matrix s->K
	*/
	// print_triplet_form("%.12g\n",t);
	if ((s->K = sparsem_from_triplet_form(t,ADDITION)) == NULL)
		goto error_malloc;

	error = 0;

error_malloc:
	if (error == 1)
		freeSystemOfEquations(s);
	freeTripletForm(t);	
	return s;
}

int writeFEM2DSolutionTXTType(DOUBLE *s,Specification2D spec)
{
	Lagrange1DValues lv1;
	Lagrange2DValues lv2;
	FILE *fp;
	DataMesh mesh;
	int tn, i, j, k, l, degree, elementvalues, ok;

	if ((lv1 = LagrangeAtPointsInBasicInterval(spec->degree,spec->elementvalues)) == NULL)
		goto final;
	if ((lv2 = LagrangeAtPointsInBasicTriangle(spec->degree,spec->elementvalues)) == NULL)
		goto final;

	if ((fp = fopen(spec->filename,"w")) == NULL)
		goto final;

	mesh = spec->mesh;
	elementvalues = spec->elementvalues;
	degree = spec->degree;
	if (elementvalues < degree + 1)
		elementvalues = degree + 1;

	/*
		Node points (elementvalues == 2)
	*/
	for(i = 0;i < mesh->nnodes;i++)
	{
		PNode nd = &(mesh->nodes[i]);
		fprintf(fp,"%.16g %.16g %.16g\n",nd->x,nd->y,s[nd->nodeno]);
	}

	if (elementvalues < 3)
		goto final;

	/*
		Interior points of the edges (elementvalues >= 3)
	*/
	for(i = 0;i < mesh->nedges;i++)
	{
		PEdge eg = &(mesh->edges[i]);
		for (j = 1; j < elementvalues - 1;j++)
		{
			node r;
			node_coordinates_at_edge(eg,elementvalues - 1,j,&r);
			DOUBLE z = 0.0;
			for(k = 0;k <= degree;k++)
			{
				int ik = global_node_number_at_edge(mesh,eg,degree,k);
				z += s[ik] * lv1[k][j].value;
			}
			fprintf(fp,"%.16g %.16g %.16g\n",r.x,r.y,z);
		}
	}

	if (elementvalues < 4)
		goto final;

	/*
		Only interior points of the triangle (elementvalues >= 4)
	*/
	for (tn = 0; tn < mesh->ntriangles; tn++)
	{
		PTriangle tr = &(mesh->triangles[tn]);
		int count = 0;
		for(k = 0;k < elementvalues;k++)
		{
			for (l = 0;l < elementvalues - k;l++)
			{
				if ((k == 0) || (l == 0) || (k + l == elementvalues - 1))
				{
					count++;
					continue;
				}
				node r;
				node_coordinates(tr,elementvalues - 1,k,l,&r);
				DOUBLE z = 0.0;
				for(i = 0;i <= degree;i++)
				{
					for(j = 0;j <= degree - i;j++)
					{
						int ij = global_node_number(mesh,tr,degree,i,j);
						z += s[ij] * lv2[i][j][count].value;
					}
				}
				fprintf(fp,"%.16g %.16g %.16g\n",r.x,r.y,z);
				count++;
			}
		}
	}

final:
	freeLagrange1DValues(lv1);
	freeLagrange2DValues(lv2,spec->degree);
	if (fp != NULL)
		fclose(fp);
	return 1;
}

int writeFEM2DPETScSolutionTXTType(PetscScalar *s,Specification2D spec)
{
	Lagrange1DValues lv1;
	Lagrange2DValues lv2;
	FILE *fp;
	DataMesh mesh;
	int tn, i, j, k, l, degree, elementvalues, ok;

	if ((lv1 = LagrangeAtPointsInBasicInterval(spec->degree,spec->elementvalues)) == NULL)
		goto final;
	if ((lv2 = LagrangeAtPointsInBasicTriangle(spec->degree,spec->elementvalues)) == NULL)
		goto final;

	if ((fp = fopen(spec->filename,"w")) == NULL)
		goto final;

	mesh = spec->mesh;
	elementvalues = spec->elementvalues;
	degree = spec->degree;
	if (elementvalues < degree + 1)
		elementvalues = degree + 1;

	/*
		Node points (elementvalues == 2)
	*/
	for(i = 0;i < mesh->nnodes;i++)
	{
		PNode nd = &(mesh->nodes[i]);
		fprintf(fp,"%.16g %.16g %.16g\n",nd->x,nd->y,s[nd->nodeno]);
	}

	if (elementvalues < 3)
		goto final;

	/*
		Interior points of the edges (elementvalues >= 3)
	*/
	for(i = 0;i < mesh->nedges;i++)
	{
		PEdge eg = &(mesh->edges[i]);
		for (j = 1; j < elementvalues - 1;j++)
		{
			node r;
			node_coordinates_at_edge(eg,elementvalues - 1,j,&r);
			DOUBLE z = 0.0;
			for(k = 0;k <= degree;k++)
			{
				int ik = global_node_number_at_edge(mesh,eg,degree,k);
				z += s[ik] * lv1[k][j].value;
			}
			fprintf(fp,"%.16g %.16g %.16g\n",r.x,r.y,z);
		}
	}

	if (elementvalues < 4)
		goto final;

	/*
		Only interior points of the triangle (elementvalues >= 4)
	*/
	for (tn = 0; tn < mesh->ntriangles; tn++)
	{
		PTriangle tr = &(mesh->triangles[tn]);
		int count = 0;
		for(k = 0;k < elementvalues;k++)
		{
			for (l = 0;l < elementvalues - k;l++)
			{
				if ((k == 0) || (l == 0) || (k + l == elementvalues - 1))
				{
					count++;
					continue;
				}
				node r;
				node_coordinates(tr,elementvalues - 1,k,l,&r);
				DOUBLE z = 0.0;
				for(i = 0;i <= degree;i++)
				{
					for(j = 0;j <= degree - i;j++)
					{
						int ij = global_node_number(mesh,tr,degree,i,j);
						z += s[ij] * lv2[i][j][count].value;
					}
				}
				fprintf(fp,"%.16g %.16g %.16g\n",r.x,r.y,z);
				count++;
			}
		}
	}

final:
	freeLagrange1DValues(lv1);
	freeLagrange2DValues(lv2,spec->degree);
	if (fp != NULL)
		fclose(fp);
	return 1;
}
