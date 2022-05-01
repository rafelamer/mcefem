#include <tfgfem.h>
#include <cmdline.h>

int main(int argc, char **argv)
{
	int ret;
	char *infile = NULL;
	char *asyfile = NULL;
	char *solver = NULL;
	int dimension;
	DOUBLE *s1 = NULL;
	PetscScalar *s2 = NULL;
	gsl_vector *s3 = NULL;
	Specification1D spec1d = NULL;
	Specification2D spec2d = NULL;
	SystemOfEquations system = NULL;

	ret = EXIT_FAILURE;

	static struct gengetopt_args_info ai;
	if (cmdline_parser(argc, argv, &ai) != 0) {
		fprintf(stderr, "Error reading the command line parameters\n");
		goto final;
	}
	if (ai.help_given) {
		printf("%s\n", gengetopt_args_info_usage);
		printf("%s\n", *gengetopt_args_info_help);
		ret = EXIT_SUCCESS;
		goto final;
	}

	infile = ai.infile_arg;
	dimension = ai.dimension_arg;
	solver = ai.solver_arg;
	if (ai.asyfile_given)
		asyfile = ai.asyfile_arg;

	if (dimension == 1)
	{
		if ((spec1d = parseXMLSpec1DDocument(infile,NULL)) == NULL)
			goto final;

		if ((system = StiffnessMatrixAndLoadVector1D(spec1d)) == NULL)
			goto final;


		// print_sparsem_matrix("%.8g  ",system->K);
		// print_vector("%.18g\n",system->F,system->K->cols);

		if (!strcmp("umfpack",solver))
		{
			if ((s1 = umfpack_solve(system->K,system->F)) == NULL)
				goto final;

			// print_vector("%.18g\n",s1,system->K->cols);

			if ( !strcmp("txt",spec1d->type))
				if ( !writeFEM1DSolutionTXTType(s1,spec1d))
					goto final;

			ret = EXIT_SUCCESS;
			goto final;
		}

		if (!strcmp("gsl",solver))
		{
			if ((s3 = gsl_gmres_solve(system)) == NULL)
				goto final;

			if ( !strcmp("txt",spec1d->type))
				if ( !writeFEM1DSolutionTXTType(gsl_vector_ptr(s3,0),spec1d))
					goto final;

			ret = EXIT_SUCCESS;
			goto final;
		}

		if (!strcmp("petsc",solver))
		{
			if (petsc_solve(system,&s2) != 0)
				goto final;

			if ( !strcmp("txt",spec1d->type))
				if ( !writeFEM1DPETScSolutionTXTType(s2,spec1d))
					goto final;
			ret = EXIT_SUCCESS;
			goto final;
		}
		goto final;
	}

	if ((spec2d = parseXMLSpec2DDocument(infile,NULL)) == NULL)
		goto final;

	if (asyfile != NULL)
		mesh_to_asy(spec2d->mesh,spec2d->degree,asyfile,1.0,0);

	if ((system = StiffnessMatrixAndLoadVector2D(spec2d)) == NULL)
		goto final;

	if (!strcmp("umfpack",solver))
	{
		if ((s1 = umfpack_solve(system->K,system->F)) == NULL)
			goto final;
		if ( !strcmp("txt",spec2d->type))
			if ( !writeFEM2DSolutionTXTType(s1,spec2d))
				goto final;
	}

	if (!strcmp("gsl",solver))
	{
		if ((s3 = gsl_gmres_solve(system)) == NULL)
			goto final;
		if ( !strcmp("txt",spec2d->type))
			if ( !writeFEM2DSolutionTXTType(gsl_vector_ptr(s3,0),spec2d))
				goto final;
	}

	if (!strcmp("petsc",solver))
	{
		if (petsc_solve(system,&s2) != 0)
			goto final;

		if ( !strcmp("txt",spec2d->type))
			if ( !writeFEM2DPETScSolutionTXTType(s2,spec2d))
				goto final;
	}

	ret = EXIT_SUCCESS;

final:
	if (ret != EXIT_SUCCESS)
		printf("Something wrong ocurred during the solution of the problem\n");
	freeSpecification1D(spec1d);
	freeSpecification2D(spec2d);
	freeSystemOfEquations(system);
	if (s1 != NULL)
		free(s1);
	if (s2 != NULL)
		PetscFree(s2);
	if (s3 != NULL)
		gsl_vector_free(s3);
	cmdline_parser_free(&ai);
	return ret;
}
