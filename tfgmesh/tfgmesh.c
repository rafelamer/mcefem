#include <tfgfem.h>
#include <cmdline.h>

int main(int argc, char **argv)
{
	int ret;
	char *infile = NULL;
	char *asyfile = NULL;
	double size;
	double unitsize;
	int degree;
	int labels;
	DataRegion rg = NULL;
	DataMesh mesh = NULL;

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
	if (ai.asyfile_given)
		asyfile = ai.asyfile_arg;
	degree = ai.degree_arg;
	size = ai.size_arg;
	unitsize = ai.unitsize_arg;
	
	if (size <= 0.0)
		goto final;

	labels = NOLABELS;
	if (ai.labels_flag)
		labels = WITHLABELS;

	if (degree < 0)
		degree = 0;

	if ((rg = parseXMLRegionDocument(infile)) == NULL) {
		printf("Error reading the file %s\n",infile);
		goto final;
	}

	size = ai.size_arg;

	if ((mesh = make_mesh(rg,size)) == NULL) {
		printf("Error generating the mesh for the file %s\n",infile);
		goto final;
	}

	if (asyfile != NULL)
		mesh_to_asy(mesh,degree,asyfile,unitsize,labels);

	if (ai.info_flag)
	{
		print_mesh_data(mesh,degree);
	}

	ret = EXIT_SUCCESS;
	
final:
	if (ret != EXIT_SUCCESS)
		printf("Something wrong ocurred during the generation of the mesh\n");

	if (mesh != NULL)
		freeDataMesh(mesh);
	if (rg != NULL)
		freeDataRegion(rg);

	cmdline_parser_free(&ai);
	return ret;
}
