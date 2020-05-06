/**************************************************************************************
* Filename:   tfgfem.h
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
#ifndef H_TFGFEM_H_
#define H_TFGFEM_H_ 1

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <matheval.h>
#include <suitesparse/umfpack.h>
#include <libxml/parser.h>
#include <libxml/xmlmemory.h>
#include <lua5.2/lua.h>
#include <lua5.2/lauxlib.h>
#include <lua5.2/lualib.h>
#include <uthash.h>
#include <Python.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_vector.h>
#include <petscksp.h>

#define DOUBLE double
#define max(a,b)            (((a) > (b)) ? (a) : (b))

/*************************************************************
 *
 *  Utility functions for working with vectors and matrices
 *
 *************************************************************/
#define FREE_APPENDED 1
#define NOT_FREE_APPENDED 0

#define make_vector(v,n) if (((v) = calloc((n),sizeof(*(v)))) == NULL)              \
	goto error_malloc;

#define free_vector(v) do { free(v); v = NULL; } while (0)

#define clone_vector(u,v,n) do {                                                    \
	make_vector(v,(n));                                                             \
	memcpy(v,(u),(n)*sizeof(*(u)));                                                 \
} while (0)

#define expand_vector(v,n) if (((v) = realloc((v),(n) * sizeof(*(v)))) == NULL)     \
	goto error_malloc;

#define append_to_vector(u,v,m,n,f) do {                                            \
	if ((u) == NULL)                                                                \
	{                                                                               \
		make_vector((u),(n));                                                       \
	}                                                                               \
	else                                                                            \
	{                                                                               \
		expand_vector((u),(m)+(n));                                                 \
	}                                                                               \
	void *p = u + m;                                                                \
	memcpy(p,(v),(n)*sizeof(*(v)));                                                 \
	if((f))                                                                         \
		free_vector((v));                                                           \
} while (0)

#define make_matrix(a, m, n) do {                                                   \
	make_vector(a, (m) + 1);                                                        \
	for (size_t make_matrix_loop_counter = 0;                                       \
		make_matrix_loop_counter < (m);                                             \
		make_matrix_loop_counter++)                                                 \
		make_vector((a)[make_matrix_loop_counter], (n));                            \
} while (0)

#define clone_matrix(u,v,m,n) do {		                                            \
	make_matrix(v,(m),(n));                                                         \
	for (size_t make_matrix_loop_counter = 0;                                       \
			make_matrix_loop_counter < (m);                                         \
			make_matrix_loop_counter++)                                             \
		memcpy(v[make_matrix_loop_counter],                                         \
			u[make_matrix_loop_counter],                                            \
			(n)*sizeof(**(v)));                                                     \
} while (0)

#define free_matrix(a) do {                                                         \
	if (a != NULL) {                                                                \
		for (size_t make_matrix_loop_counter = 0;                                   \
				(a)[make_matrix_loop_counter] != NULL;                              \
				make_matrix_loop_counter++)                                         \
			free_vector((a)[make_matrix_loop_counter]);                             \
		free_vector(a);                                                             \
		a = NULL;                                                                   \
	}                                                                               \
} while (0)

#define print_vector(fmt, v, n) do {                                                \
	for (size_t print_vector_loop_counter = 0;                                      \
			print_vector_loop_counter < (n);                                        \
			print_vector_loop_counter++)                                            \
		printf(fmt, (v)[print_vector_loop_counter]);                                \
	putchar('\n');                                                                  \
} while (0)

#define print_matrix(fmt, a, m, n) do {                                             \
	for (size_t print_matrix_loop_counter = 0;                                      \
			print_matrix_loop_counter < (m);                                        \
			print_matrix_loop_counter++)                                            \
		print_vector(fmt, (a)[print_matrix_loop_counter], (n));                     \
} while (0)

#define make_triangular_matrix(a, n) do {                                           \
	make_vector(a, (n) + 1);                                                        \
	for (size_t make_matrix_loop_counter = 0;                                       \
			make_matrix_loop_counter < (n);                                         \
			make_matrix_loop_counter++)                                             \
		make_vector((a)[make_matrix_loop_counter],                                  \
			(n) - make_matrix_loop_counter );		                                \
	(a)[n] = NULL;                                                                  \
} while (0)  

#define matrix_to_vector(a,v,m,n,type) do {                                         \
	type *matrix_vector_pointer;                                                    \
	matrix_vector_pointer = v;			                                            \
	for (size_t matrix_vector_loop_counter = 0;                                     \
			matrix_vector_loop_counter < (n);                                       \
			matrix_vector_loop_counter++) {		                                    \
		memcpy(matrix_vector_pointer,                                               \
		a[matrix_vector_loop_counter],(n)*sizeof(type));                            \
		matrix_vector_pointer += (n);                                               \
	}	           		                                                            \
} while (0)

#define vector_to_matrix(v,a,m,n,type) do {                                         \
	type *vector_matrix_pointer;                                                    \
	vector_matrix_pointer = v;			                                            \
	for (size_t vector_matrix_loop_counter = 0;                                     \
			vector_matrix_loop_counter < (m);                                       \
			vector_matrix_loop_counter++) {		                                    \
		memcpy(a[vector_matrix_loop_counter],                                       \
			vector_matrix_pointer,(n)*sizeof(type));                                \
		vector_matrix_pointer += (n);                                               \
	}	           		                                                            \
} while (0)

/*************************************************************
 *
 * Read data attributes from an XML file
 *
 *************************************************************/
#define XML_VALUE_REQUIRED 1
#define XML_VALUE_OPTIONAL 0
int get_int_value(xmlNodePtr node,char *name,int *found,int optional);
double get_double_value(xmlNodePtr node,char *name,int *found,int optional);
xmlChar *get_string_value(xmlNodePtr node,char *name,int *found,int optional);
xmlChar *get_string_content(xmlDocPtr doc,xmlNodePtr node);

/*************************************************************
 * 
 * A region or domain in dimension 2 is a collections a closed
 * curves. For each curve we have to say if it represents a hole
 * in the doman.
 * Every close curve is defined by a collection of points
 *
 *************************************************************/
typedef struct {
	int point_no;
	DOUBLE x;
	DOUBLE y;
	int bc;          /* Boundary condition of the point */
	int sbc;         /* Boundary condition of the segment until the next point */
} data_point;

typedef struct {
	int segment_no;
	int point_no_1;
	int point_no_2;
	int bc;          /* Boundary condition */
} data_segment;

typedef struct {
	DOUBLE x;
	DOUBLE y;
} data_hole;

typedef struct {
	data_point *points;
	data_segment *segments;
	data_hole *holes;
	int npoints;
	int nsegments; 
	int nholes;
} data_region;
typedef data_region *DataRegion;

void free_data_region(DataRegion *r);
#define freeDataRegion(n) free_data_region(&(n))

DataRegion data_region_new();
DataRegion parseXMLRegionDocument(char *docname);
void region_to_asy(DataRegion r,const char *outfile);
void print_region(DataRegion r);

/*************************************************************
 *
 * Triangulating a planar polygonal domain means dividing the domain 
 * into a union of triangles so that
 *    1. every triangle has a nonempty interior;
 *    2. no two triangles have common interior parts; and
 *    3. any side of any triangle is either a part of the domain’s 
 *       boundary or is the side of another triangle.
 *
 * A triangulated domain is also called a triangular mesh.
 *
 *************************************************************/
#define FEM_BC_NOT_DEFINED 0
#define FEM_BC_DIRICHLET 2
#define FEM_BC_NEUMANN 3
#define FEM_BC_ROBIN 4

typedef struct {
	int nodeno;                   /* Numeration of all nodes */
	int bnodeno;                  /* Numeration of the boundary nodes */
	DOUBLE x;
	DOUBLE y;
	int bc;
} node;
typedef node *PNode;

typedef struct {
	int edgeno;                   /* Numeration of all edges */
	int bedgeno;                  /* Numeration of the boundary edges */
	node *n[2];
	int bc;
} edge;
typedef edge *PEdge;

typedef struct {
	int triangleno;
	node *n[3];
	edge *e[3];
	DOUBLE ex[3], ey[3];          /* x and y components of the edge vectors */
	DOUBLE area;                  /* Triangle’s area */
} triangle_data;
typedef triangle_data *PTriangle;

typedef struct {
	node *nodes;                  /* The array of node structures */
	edge *edges;                  /* The array of edge structures */
	triangle_data *triangles;     /* The array of triangle structures */
	int nnodes, nbnodes;
	int nedges, nbedges;
	int ntriangles;
	DOUBLE area;                  /* Total area of the region */
	int hasDirichletBC, hasNeumannBC, hasRobinBC;
} data_mesh;
typedef data_mesh *DataMesh;

void free_data_mesh(DataMesh *m);
#define freeDataMesh(n) free_data_mesh(&(n))

#define NOLABELS 0
#define WITHLABELS 1

DataMesh mesh_new();
DataMesh make_mesh(DataRegion r,DOUBLE a);
void mesh_to_asy(DataMesh mesh,int degree,const char *outfile,double unitsize,int labels);
void mesh_to_asy_boundary(DataMesh mesh,const char *outfile);
void mesh_to_asy_region(DataMesh mesh,const char *outfile);
DataMesh parseXMLMeshDocument(char *docname);
int writeXMLMeshDocument(DataMesh mesh,char *docname);
int number_of_u_unknowns(DataMesh mesh,int degree);
int number_of_dotu_unknowns(DataMesh mesh,int degree);
int number_of_unknowns(DataMesh mesh,int degree);
void print_mesh_data(DataMesh mesh,int degree);
void print_triangle(DataMesh mesh,int n);
void print_list_of_distinct_points(DataMesh mesh,int degree,const char *outfile);
void print_list_of_points_per_triangle(DataMesh mesh,int degree,const char *outfile);
int global_node_number(DataMesh mesh,PTriangle tr,int degree,int m,int n);
int global_node_number_at_edge(DataMesh mesh,PEdge edge,int degree,int m);
int boundary_node_number(DataMesh mesh,PTriangle tr,int degree,int m,int n);
int boundary_node_number_at_edge(DataMesh mesh,PEdge edge,int degree,int m);
int boundary_condition_at_node_edge(PEdge edge,int degree,int m);
void node_coordinates(PTriangle tr,int degree,int m,int n,PNode r);
void node_coordinates_at_edge(PEdge edge,int degree,int m,PNode r);

/*************************************************************
 *
 * Gauss quadrature in dimensions 1 and 2
 *
 * The points and weights are defined in the interval [-1,1] in
 * dimension 1 or in the triangle with vertices (0,0) 
 * (1,0) and (0,1) when the dimension is 2.
 *
 *************************************************************/
typedef DOUBLE (*FunctionX)(DOUBLE);
typedef DOUBLE (*FunctionXY)(DOUBLE,DOUBLE);

typedef struct {
	DOUBLE x;	/* point */
 	DOUBLE w;	/* weight */
} gauss_quadrature;
typedef gauss_quadrature *GaussQuadrature;

typedef struct {
  DOUBLE lambda1;
  DOUBLE lambda2;
  DOUBLE weight;
} twb_gauss_quadrature;
typedef twb_gauss_quadrature *TWBGaussQuadrature;

GaussQuadrature gauss_data(int *n);
DOUBLE gauss_integrate_function(DOUBLE a,DOUBLE b, FunctionX f,int points);
DOUBLE gauss_integrate_evaluator(DOUBLE a,DOUBLE b,void *evaluator,int points);

TWBGaussQuadrature get_twb_data(int *d, int *n);
DOUBLE twb_integrate_triangle_function(PTriangle tr,FunctionXY f,TWBGaussQuadrature qdat);
DOUBLE twb_integrate_triangle_evaluator(PTriangle tr,void *evaluator,TWBGaussQuadrature qdat);
DOUBLE twb_integrate_mesh_function(DataMesh mesh,FunctionXY f,int degree);
DOUBLE twb_integrate_mesh_evaluator(DataMesh mesh,void *evaluator,int degree);

/*************************************************************
 *
 *	Lagrange polynomials in dimensions 1 and 2
 *  They are defined in the interval [0,1] or in the triangle
 *  with vertices (0,0) (1,0) and (0,1)
 * 
 *************************************************************/
typedef struct {
	DOUBLE value; 
	DOUBLE dvalue; /* derivative of Lagrange polynomial */
} lagrange_1d_values;
typedef lagrange_1d_values **Lagrange1DValues;

typedef struct {
	DOUBLE value;
	DOUBLE dxvalue;
	DOUBLE dyvalue;
} lagrange_2d_values;
typedef lagrange_2d_values ***Lagrange2DValues;

void free_lagrange_1d_values(Lagrange1DValues *data);
#define freeLagrange1DValues(n) free_lagrange_1d_values(&(n))
void free_lagrange_2d_values(Lagrange2DValues *data,int elements);
#define freeLagrange2DValues(n,d) free_lagrange_2d_values(&(n),(d))

DOUBLE Lagrange1D(DOUBLE x,unsigned int i,unsigned int n);
DOUBLE DerivativeLagrange1D(DOUBLE x,unsigned int i,unsigned int n);
DOUBLE Lagrange2D(DOUBLE x,DOUBLE y,unsigned int i,unsigned j,unsigned int n);
DOUBLE DerivativeXLagrange2D(DOUBLE x,DOUBLE y,unsigned int i,unsigned j,unsigned int n);
DOUBLE DerivativeYLagrange2D(DOUBLE x,DOUBLE y,unsigned int i,unsigned j,unsigned int n);

Lagrange1DValues LagrangeAtGaussPoints(int degree,GaussQuadrature qdat,int npoints);
Lagrange2DValues LagrangeAtTWBPoints(int degree,TWBGaussQuadrature qdat,int npoints);

Lagrange1DValues LagrangeAtPointsInBasicInterval(int degree,int npoints);
Lagrange2DValues LagrangeAtPointsInBasicTriangle(int degree,int npoints);

/*************************************************************
 *
 * Sparse matrices
 * Programming Projects in C for Students of Engineering, Science, and Mathematics
 * Rouben Rostamian
 * ISBN 978-1-611973-49-5
 * SIAM 2015
 *
 * Direct Methods for Sparse Linear Systems
 * Timothy A. Davis
 * ISBN: 978-0898716139
 * SIAM 2006 
 *
 *************************************************************/
typedef struct {
  unsigned int rows;       /* Number of rows                                        */
  unsigned int cols;       /* Number of columns                                     */
  unsigned int nz;         /* Number of non zero entries                            */
  unsigned int size_of_Ax; /* Size allocated for vectors Ax and Ai                  */
  DOUBLE *Ax;              /* Non zero entries of the matrix                        */
  unsigned int *Ap;        /* Number of nonzero elements in the first i columns     */
  unsigned int *Ai;        /* Row index in the orginal matrix of each element of Ax */
} sparse_matrix;
typedef sparse_matrix *SparseMatrix;

typedef struct { 
  unsigned int row;
  unsigned int col;
  DOUBLE value;
} MatrixElement;

typedef struct {
  unsigned int elements;
  unsigned int size_of_T;
  unsigned int rows;
  unsigned int cols;
  MatrixElement *T;
} triplet_form;
typedef triplet_form *TripletForm;

#define freeSparseMatrix(n)   free_sparsem_matrix(&(n))
#define freeTripletForm(n)    free_triplet_form(&(n))
#define SIZEINCREMENT 1024
#define REPLACE 0
#define ADDITION 1

DOUBLE **sparsem_unpack(SparseMatrix r);
SparseMatrix sparsem_init(unsigned int rows,unsigned int cols,unsigned int size);
SparseMatrix sparsem_pack(DOUBLE **m,unsigned int rows,unsigned int cols);
DOUBLE matrix_element(SparseMatrix r,unsigned int row,unsigned int col);
void print_sparsem_matrix(const char *fmt,SparseMatrix r);
void print_sparse_matrix_and_vector_to_txt_file(SparseMatrix r,DOUBLE *v,unsigned int length,const char *filename);
void print_sparsem_matrix_elements(const char *fmt,SparseMatrix r);

void sparsem_write(FILE *fp,SparseMatrix r);
SparseMatrix sparsem_read(FILE *fp);

int sparsem_change_element(SparseMatrix r,DOUBLE value,unsigned int row,unsigned int col,int mode);
DOUBLE *sparsem_multiply_by_vector(SparseMatrix r,DOUBLE *x);
SparseMatrix sparsem_from_triplet_form(TripletForm t,int mode);
TripletForm triplet_form_init(unsigned int size);
TripletForm triplet_form_read(const char *filename,unsigned char symmetric,unsigned int zero_based);
int triplet_form_append_element(TripletForm t,DOUBLE value,unsigned int row,unsigned int col);
void print_triplet_form(const char *fmt,TripletForm t);

void free_sparsem_matrix(SparseMatrix *r);
void free_triplet_form(TripletForm *t);
DOUBLE *umfpack_solve(SparseMatrix r,DOUBLE *b);

/*************************************************************
 *
 * Different types of functions used. 
 *   1. Matheval evaluators
 *   2. C functions
 *   3. Lua functions
     4. Pyton functions
 *
**************************************************************/
#define MATHEVALUATOR 1
#define CFUNCTION1D 2
#define LUAFUNCTION 3
#define CFUNCTION2D 4
#define PYTHONFUNCTION 5

typedef struct {
	char name[128];
	int type;
	lua_State *Lua;
	PyObject *Python;
	void *evaluator;
	FunctionX f1d;
	FunctionXY f2d;
	UT_hash_handle hh;
} function_table;
typedef function_table *FunctionTable;

PyObject *Py_Initialize_Functions(char *programname,char *filename);
void Py_Finalize_Functions(PyObject *module);

void add_pythonfunction(FunctionTable *hash,char *name,PyObject *python);
void add_cfunction1d(FunctionTable *hash,char *name,FunctionX f);
void add_cfunction2d(FunctionTable *hash,char *name,FunctionXY f);
void add_mathevaluator(FunctionTable *hash,char *name,char *function);
void add_luafunction(FunctionTable *hash,char *name,lua_State *Lua,char *function);
DOUBLE callFunctionTable1D(FunctionTable hash,char *name,DOUBLE x);
DOUBLE callFunctionTable2D(FunctionTable hash,char *name,DOUBLE x,DOUBLE y);
DOUBLE runFunctionTable1D(FunctionTable item,DOUBLE x);
DOUBLE runFunctionTable2D(FunctionTable item,DOUBLE x,DOUBLE y);
void free_function_table(FunctionTable *table);
#define freeFunctionTable(n)  free_function_table(&(n))

/*************************************************************
 *
 * Finite Element Method in dimension 1. 
 * Stiffness matrix K and load vector F
 *
 *************************************************************/
typedef struct {
	DOUBLE xa, xb;                /* Interval [xa,xb] */
	int elements;                 /* Number of subintervals or elements */
	int degree;                   /* Degree of the Lagrange Polynomials */
	Lagrange1DValues lv;          /* Lagrage values at Gauss points */
	int npoints;                  /* Number of points for Gauss quadrature */
	GaussQuadrature qdat;         /* Points and weights for Gauss quadrature */ 
	FunctionTable funs;           /* Table for functions a2, a1, a0 and f */
	int bctype[2];                /* Type of boundary conditions at xa and xb */
	DOUBLE bc[2][3];              /* Boundary conditions */
	char type[64];                /* Type of output */
	char filename[1024];          /* Output file */
	PyObject *pModule;            /* Python object module */
	lua_State *Lua;               /* Lua interpreter */
	int elementvalues;            /* Output values per element */
} spec1D;
typedef spec1D *Specification1D;

void free_spec1D(Specification1D *spec);
#define freeSpecification1D(n)   free_spec1D(&(n))

Specification1D parseXMLSpec1DDocument(char *docname,FunctionTable funs);

typedef struct {
	SparseMatrix K;
	DOUBLE *F;
} system_of_equations;
typedef system_of_equations *SystemOfEquations;

void free_system_of_equations(SystemOfEquations *s);
#define freeSystemOfEquations(n)   free_system_of_equations(&(n))

SystemOfEquations StiffnessMatrixAndLoadVector1D(Specification1D spec);
int writeFEM1DSolutionTXTType(DOUBLE *s,Specification1D spec);
int writeFEM1DPETScSolutionTXTType(PetscScalar *s,Specification1D spec);

/*************************************************************
 *
 * Finite Element Method in dimension 2. 
 * Stiffness matrix K and load vector F
 *
 *************************************************************/
#define NODERIVATIVE 0
#define DXDERIVATIVE 1
#define DYDERIVATIVE 2
#define NODERIVATIVES 3
#define DXNODERIVATIVES 4
#define NODXDERIVATIVES 5
#define DYNODERIVATIVES 6
#define NODYDERIVATIVES 7
#define DXDXDERIVATIVES 8
#define DYDYDERIVATIVES 9
#define DXDYDERIVATIVES 10
#define DYDXDERIVATIVES 11
#define NABLADERIVATIVES 12

DOUBLE integralOverOneTriangleWithOneLagrange(PTriangle tr,unsigned int m1,unsigned int n1,unsigned int d1,
											  FunctionTable item,TWBGaussQuadrature q2,Lagrange2DValues lv);

DOUBLE integralOverOneTriangleWithTwoLagrange(PTriangle tr,unsigned int m1,unsigned int n1,
											  unsigned int m2,unsigned int n2,unsigned int derivative,  
											  FunctionTable item,TWBGaussQuadrature q2,Lagrange2DValues lv);

DOUBLE integralOverOneTriangleWithTwoLagrangeTest(PTriangle tr,unsigned int m1,unsigned int n1,
												  unsigned int m2,unsigned int n2,
												  FunctionTable item,TWBGaussQuadrature q2,Lagrange2DValues lv2);

DOUBLE integralFunctionOverOneEdge(PTriangle tr,int edge,FunctionTable item,GaussQuadrature q1);

DOUBLE integralOverOneEdgeTwoLagrange(PTriangle tr,int edge,unsigned int m1,unsigned int n1,
									  FunctionTable item,GaussQuadrature q1,Lagrange1DValues lv);

typedef struct {
	DataMesh mesh;                /* The triangulation and boundary elements */
	int degree;                   /* Degree of the Lagrange Polynomials */
	int npoints1d;                /* Number of points for Gauss quadrature in 1D */
	Lagrange1DValues lv1d;        /* Lagrage values at Gauss points */
	GaussQuadrature qdat1d;       /* Points and weights for Gauss quadrature in 1D */
	int npoints2d;                /* Number of points for TWB quadrature in 2D */
	Lagrange2DValues lv2d;        /* Lagrage values at TWB points */
	TWBGaussQuadrature qdat2d;    /* Points and weights for TWB quadrature in 2D */
	FunctionTable funs;           /* Table for functions a, b1, b2, c, f, d, n, A, B and C */
	char type[64];                /* Type of output */
	char filename[1024];          /* Output file */
	PyObject *pModule;            /* Python object module */
	lua_State *Lua;               /* Lua interpreter */
	int elementvalues;            /* Output values per element */
} spec2D;
typedef spec2D *Specification2D;

void free_spec2D(Specification2D *spec);
#define freeSpecification2D(n)   free_spec2D(&(n))

Specification2D parseXMLSpec2DDocument(char *docname,FunctionTable funs);
SystemOfEquations StiffnessMatrixAndLoadVector2D(Specification2D spec);

int writeFEM2DSolutionTXTType(DOUBLE *s,Specification2D spec);
int writeFEM2DPETScSolutionTXTType(PetscScalar *s,Specification2D spec);

/*************************************************************
 *
 * Iterative solver of a system linear equations by the
 * Generalized Minimum Residual Method (GMRES) using
 * the GNU GSL library
 *
 *************************************************************/

gsl_vector *gsl_gmres_solve(SystemOfEquations system);

/*************************************************************
 *
 * PETSC solver
 *
 *************************************************************/

PetscScalar *petsc_solve(SystemOfEquations system);

#endif
