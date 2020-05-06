/**************************************************************************************
* Filename:   gsl.c
* Authors:     
* Copyright:  
* Disclaimer: This code is presented "as is" and it has been written to 
*             implement the Finite Element Method in dimensions 1 and 2.
*             It has been writen for educational purposes.
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
#include <gsl/gsl_splinalg.h>

gsl_vector *gsl_gmres_solve(SystemOfEquations system)
{
	unsigned int size;
	size_t i, j;


	size = system->K->cols;
	gsl_spmatrix *A = gsl_spmatrix_alloc(size,size);      /* triplet format */
	gsl_spmatrix *C;                                      /* compressed format */
	gsl_vector *f = gsl_vector_alloc(size);               /* right hand side vector */
	gsl_vector *u = gsl_vector_alloc(size);               /* solution vector */

	/*
		Set the rigth hand side vector
	*/
	for (i = 0; i < size;i++)
		gsl_vector_set(f, i, system->F[i]);

	/*
		Set the sparse matrix
	*/
	for (i = 0;i < size;i++)
		for (j = system->K->Ap[i];j < system->K->Ap[i+1];j++)
			gsl_spmatrix_set(A,system->K->Ai[j],i,system->K->Ax[j]);

	C = gsl_spmatrix_ccs(A);
	const double tol = 1.0e-15;                                /* solution relative tolerance */
    const size_t max_iter = 10000000;                          /* maximum iterations */
	const gsl_splinalg_itersolve_type *T = gsl_splinalg_itersolve_gmres;
	gsl_splinalg_itersolve *work = gsl_splinalg_itersolve_alloc(T,size,0);

	size_t iter = 0;
    double residual;
    int status;

	gsl_vector_set_zero(u);
	do 
		status = gsl_splinalg_itersolve_iterate(C,f,tol,u,work);
    while (status == GSL_CONTINUE && ++iter < max_iter);

	gsl_splinalg_itersolve_free(work);
	gsl_spmatrix_free(A);
	gsl_spmatrix_free(C);
	gsl_vector_free(f);

    if (status != GSL_SUCCESS)
    {
    	gsl_vector_free(u);
    	printf("Impossible to solve the system of equations after %lu iterations\n",iter);
    	u = NULL;
    }
    else
    {
    	printf("Solved the system of equations after %lu iterations\n",iter);
    }
	return u;
}