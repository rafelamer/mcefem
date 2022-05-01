/**************************************************************************************
* Filename:   petsc.c
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

PetscErrorCode petsc_solve(SystemOfEquations system,PetscScalar **s)
{
	Vec x, B;
	Mat A;
	KSP ksp;
	PC  pc;
	PetscErrorCode ierr;
	PetscInt i, j, size;
	PetscInt *idx;
	static char help[] = "Solves a linear system with PETSc\n\n";

	PetscInitialize(NULL,NULL,(char*)0,help);
	size = system->K->rows;
	ierr = PetscMalloc1(size,s); CHKERRQ(ierr);
	ierr = PetscMalloc1(size,&idx); CHKERRQ(ierr);
	for (i = 0;i < size;i++)
    	idx[i] = i;

  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject) x, "Solution");CHKERRQ(ierr);
  ierr = VecSetSizes(x,PETSC_DECIDE,size);CHKERRQ(ierr);
	ierr = VecSetFromOptions(x);CHKERRQ(ierr);
	ierr = VecDuplicate(x,&B);CHKERRQ(ierr);

	/*
		Copy system->F to B
	*/
	ierr = VecSetValues(B,size,idx,system->F,INSERT_VALUES);CHKERRQ(ierr);

	/*
		Copy system->K to A
	*/
	ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
	ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,size,size);CHKERRQ(ierr);
	ierr = MatSetFromOptions(A);CHKERRQ(ierr);
	ierr = MatSetUp(A);CHKERRQ(ierr);

	for (i = 0;i < size;i++)
		for (j = system->K->Ap[i];j < system->K->Ap[i+1];j++)
		{
			ierr = MatSetValue(A,system->K->Ai[j],i,system->K->Ax[j],INSERT_VALUES);
			CHKERRQ(ierr);
		}
	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);

	ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
	ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);
	ierr = KSPSetTolerances(ksp,1.e-25,PETSC_DEFAULT,PETSC_DEFAULT,1000000);CHKERRQ(ierr);

	ierr = KSPSolve(ksp,B,x);CHKERRQ(ierr);

	ierr = VecGetValues(x,size,idx,*s);CHKERRQ(ierr);

	ierr = VecDestroy(&x);CHKERRQ(ierr);
	ierr = VecDestroy(&B);CHKERRQ(ierr);
	ierr = MatDestroy(&A);CHKERRQ(ierr);
	ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
	ierr = PetscFree(idx);CHKERRQ(ierr);

	ierr = PetscFinalize();

	return 0;
}
