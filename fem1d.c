/**************************************************************************************
* Filename:   fem1d.c
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
#include <string.h>
#include <ctype.h>

void free_system_of_equations(SystemOfEquations *s)
{
	if (*s == NULL)
		return;
	freeSparseMatrix((*s)->K);
	if ((*s)->F != NULL) 
		free_vector((*s)->F);
	free(*s);
	*s = NULL;
}

SystemOfEquations StiffnessMatrixAndLoadVector1D(Specification1D spec)
{
	SystemOfEquations s = NULL;
	TripletForm t = NULL;
	Lagrange1DValues lv;
	int error = 1;
	int unknowns, k, i, j, degree;
	FunctionTable a2, a1, a0, f;

	/*
		The number of unknowns is
		The number of points: spec->elements * spec->degree + 1 
		plus 2: u'(spec->xa) and u'(spec->xb)  
	*/
	unknowns = spec->elements * spec->degree + 3;
	degree = spec->degree;
	lv = spec->lv;

	HASH_FIND_STR(spec->funs,"a2",a2);
	HASH_FIND_STR(spec->funs,"a1",a1);
	HASH_FIND_STR(spec->funs,"a0",a0);
	HASH_FIND_STR(spec->funs,"f",f);

	if ((a2 == NULL) || (a1 == NULL) || (a0 == NULL) || (f == NULL)) 
		goto error_malloc;

	if ((s = (SystemOfEquations)malloc(sizeof(system_of_equations))) == NULL)
		goto error_malloc;

	if ((t = triplet_form_init(SIZEINCREMENT)) == NULL)
		goto error_malloc;

	make_vector(s->F,unknowns);

	/*
		Generation of t and s->F
	*/
	DOUBLE h = (spec->xb - spec->xa) / spec->elements;
	for (k = 0;k < spec->elements;k++)
	{
		DOUBLE xa, xb, alpha, beta;
		xa = spec->xa + k * h;
		xb = xa + h;
		alpha = 0.5 * (xb - xa);
		beta = 0.5 * (xb + xa);
		for (i = 0; i <= degree;i++)
		{
			DOUBLE r, x;
			GaussQuadrature q;
			int n;
			/*
				Computation of the load vector s->F[degree * k + i]
				r1 = \int_{xa}^{xb} f(x) \phi_{degree * k + i}(x) dx 
			*/
			q = spec->qdat;
			r = 0.0;
			n = 0;
			while(q->w != -1)
			{
				x = alpha * q->x + beta;
				r += q->w * runFunctionTable1D(f,x) * lv[i][n].value;
				n++;
				q++;
			}
			s->F[degree * k + i] += alpha * r;

			for (j = 0;j <= i;j++)
			{
				// printf("Computing element (%d,%d)\n",degree * k + i,degree * k + j);
				/*
					Computation of the stiffness matrix triplet form

					r1 =  \int_{xa}^{xb} a_0(x) \phi_{degree*k+i}(x) \phi_{degree*k+j}(x) dx    Symmetric respect i and j

					r2 = \int_{xa}^{xb} a_2(x) \phi'_{degree*k+i}(x) \phi'_{degree*k+j}(x) dx    Symmetric respect i and j

					r3 = \int_{xa}^{xb} a_1(x) \phi_{degree*k+i}(x) \phi'_{degree*k+j}(x) dx

					r4 = \int_{xa}^{xb} a_1(x) \phi'_{degree*k+i}(x) \phi_{degree*k+j}(x) dx
				*/
				DOUBLE r1, r2, r3, r4;
				n = 0;
				q = spec->qdat;
				r1 = r2 = r3 = r4 = 0.0;
				while(q->w != -1)
				{
					x = alpha * q->x + beta;
					r1 += q->w * runFunctionTable1D(a0,x) * lv[i][n].value * lv[j][n].value;
					r2 += q->w * runFunctionTable1D(a2,x) * lv[i][n].dvalue * lv[j][n].dvalue;
					r3 += q->w * runFunctionTable1D(a1,x) * lv[i][n].value * lv[j][n].dvalue;
					if (i != j)
						r4 += q->w * runFunctionTable1D(a1,x) * lv[i][n].dvalue * lv[j][n].value;
					n++;
					q++;
				}
				r1 *= alpha;
				r2 /= (4.0 * alpha);
				r3 /= 2.0;
				r4 /= 2.0;
				if (i == j)
				{
					if (! triplet_form_append_element(t,r1+r2,degree * k + i,degree * k + i))
						goto error_malloc;
					if (r3 != 0)
						if (! triplet_form_append_element(t,r3,degree * k + i,degree * k + i))
							goto error_malloc;
				}
				else
				{
					if (! triplet_form_append_element(t,r1+r2,degree * k + i,degree * k + j))
						goto error_malloc;
					if (! triplet_form_append_element(t,r1+r2,degree * k + j,degree * k + i))
						goto error_malloc;
					if (r3 != 0.0)
						if (! triplet_form_append_element(t,r3,degree * k + i,degree * k + j))
							goto error_malloc;
					if (r4 != 0.0)
						if (! triplet_form_append_element(t,r4,degree * k + j,degree * k + i))
							goto error_malloc;
				}
			} 
		}
	}
	
	/*
		Additional unknowns u'(spec->xa) and u'(spec->xb)
	 */
	if (! triplet_form_append_element(t,runFunctionTable1D(a2,spec->xa),0,unknowns - 2))
		goto error_malloc;
	if (! triplet_form_append_element(t,- runFunctionTable1D(a2,spec->xb),unknowns - 3,unknowns - 1))
		goto error_malloc;
	
	/*
		Left boundary conditions: Penultimate equation
	*/
	switch (spec->bctype[0])
	{
		case FEM_BC_DIRICHLET:
			if (! triplet_form_append_element(t,1.0,unknowns - 2,0))
				goto error_malloc;
			s->F[unknowns - 2] = spec->bc[0][0];
			break;
		case FEM_BC_NEUMANN:
			if (! triplet_form_append_element(t,1.0,unknowns - 2,unknowns - 2))
				goto error_malloc;
			s->F[unknowns - 2] = spec->bc[0][0];
			break;
		case FEM_BC_ROBIN:
			if (! triplet_form_append_element(t,spec->bc[0][0],unknowns - 2,0))
				goto error_malloc;
			if (! triplet_form_append_element(t,spec->bc[0][1],unknowns - 2,unknowns - 2))
				return 0;
			s->F[unknowns - 2] = spec->bc[0][2];
	}

	/*
		Right boundary conditions: Last equation
	*/
	switch (spec->bctype[1])
	{
		case FEM_BC_DIRICHLET:
			if (! triplet_form_append_element(t,1.0,unknowns - 1,unknowns - 3))
				goto error_malloc;
			s->F[unknowns - 1] = spec->bc[1][0];
			break;
		case FEM_BC_NEUMANN:
			if (! triplet_form_append_element(t,1.0,unknowns - 1,unknowns - 1))
				goto error_malloc;
			s->F[unknowns - 1] = spec->bc[1][0];
			break;

		case FEM_BC_ROBIN:
			if (! triplet_form_append_element(t,spec->bc[1][0],unknowns - 1,0))
				goto error_malloc;
			if (! triplet_form_append_element(t,spec->bc[1][1],unknowns - 1,unknowns - 2))
				goto  error_malloc;
			s->F[unknowns - 1] = spec->bc[1][2];
	}

	/*
		Generate the sparse matrix s->K
	*/
	if ((s->K = sparsem_from_triplet_form(t,ADDITION)) == NULL)
		goto error_malloc;

	// print_triplet_form("%.18g\n",t);
	// print_sparsem_matrix("%.8g  ",s->K);
	// print_vector("%.18g\n",s->F,unknowns);

	error = 0;

error_malloc:
	if (error == 1)
		freeSystemOfEquations(s);
	freeTripletForm(t);	
	return s; 
}

int writeFEM1DSolutionTXTType(DOUBLE *s,Specification1D spec)
{
	Lagrange1DValues lv;
	int k, i, j, ret;
	DOUBLE p, q, xa;
	FILE *fp;

	p = (spec->xb - spec->xa) / spec->elements;
	q = p / (spec->elementvalues - 1);
	ret = 0;
	fp = NULL;

	if ((lv = LagrangeAtPointsInBasicInterval(spec->degree,spec->elementvalues)) == NULL)
		goto final;

	if ((fp = fopen(spec->filename,"w")) == NULL)
		goto final;

	for (k = 0;k < spec->elements;k++)
	{
		DOUBLE x;
		xa = spec->xa + k * p;
		for (j = 0;j < spec->elementvalues - 1;j++)
		{   
			DOUBLE x, y;
			y = 0.0;
			x = xa + j * q;
			for(i = 0;i <= spec->degree;i++)
				y += s[i] * lv[i][j].value;
			if (fprintf(fp,"%.16g\t%.16g\n",x,y) < 0)
			{
				fprintf(stderr,"Error writing to file %s\n",spec->filename);
				exit(EXIT_FAILURE);
			}		
		}
		s += spec->degree; 
	}
	if (fprintf(fp,"%.16g\t%.16g\n",spec->xb,s[0]) < 0)
	{
		fprintf(stderr,"Error writing to file %s\n",spec->filename);
		exit(EXIT_FAILURE);
	}	
	ret = 1;

final:
	freeLagrange1DValues(lv);
	if (fp != NULL)
		fclose(fp);
	return ret;
}

int writeFEM1DPETScSolutionTXTType(PetscScalar *s,Specification1D spec)
{
	Lagrange1DValues lv;
	int k, i, j, ret;
	DOUBLE p, q, xa;
	FILE *fp;

	p = (spec->xb - spec->xa) / spec->elements;
	q = p / (spec->elementvalues - 1);
	ret = 0;
	fp = NULL;

	if ((lv = LagrangeAtPointsInBasicInterval(spec->degree,spec->elementvalues)) == NULL)
		goto final;

	if ((fp = fopen(spec->filename,"w")) == NULL)
		goto final;

	for (k = 0;k < spec->elements;k++)
	{
		DOUBLE x;
		xa = spec->xa + k * p;
		for (j = 0;j < spec->elementvalues - 1;j++)
		{   
			DOUBLE x, y;
			y = 0.0;
			x = xa + j * q;
			for(i = 0;i <= spec->degree;i++)
				y += s[i] * lv[i][j].value;
			if (fprintf(fp,"%.16g\t%.16g\n",x,y) < 0)
			{
				fprintf(stderr,"Error writing to file %s\n",spec->filename);
				exit(EXIT_FAILURE);
			}		
		}
		s += spec->degree; 
	}
	if (fprintf(fp,"%.16g\t%.16g\n",spec->xb,s[0]) < 0)
	{
		fprintf(stderr,"Error writing to file %s\n",spec->filename);
		exit(EXIT_FAILURE);
	}	
	ret = 1;

final:
	freeLagrange1DValues(lv);
	if (fp != NULL)
		fclose(fp);
	return ret;
}
