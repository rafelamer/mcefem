/**************************************************************************************
* Filename:   sparsem.c
* Author:     
* Copyright:  
* Disclaimer: This code is presented "as is" and it has been written to 
*             implement the Finite Element Method in dimension 2. It
*             has been writen educational purposes.
*	    
* License:    This library is free software; you can redistribute it and/or
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
*	            See https://www.gnu.org/licenses/
***************************************************************************************/
#include <tfgfem.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <suitesparse/umfpack.h>

DOUBLE **sparsem_unpack(SparseMatrix r)
{
  DOUBLE **m;
  unsigned int i, j, k;

	m = NULL;
	make_matrix(m,r->rows,r->cols);
  for (j = 0;j < r->cols;j++)
    for (k = r->Ap[j];k < r->Ap[j+1];k++)
    {
      i = r->Ai[k];
      m[i][j] = r->Ax[k];
    }
  return m;

error_malloc:
	if (m == NULL)
		return NULL;
	for (k = 0;k < r->rows;k++)
		free_vector(m[k]);
	free(m);
	return NULL;
}

SparseMatrix sparsem_init(unsigned int rows,unsigned int cols,unsigned int size)
{
  SparseMatrix r;
  if ((r = (SparseMatrix)malloc(sizeof(sparse_matrix))) == NULL)
    return NULL;
  r->rows = rows;
  r->cols = cols;
  r->nz = 0;
  r->size_of_Ax = size;
  make_vector(r->Ap,cols+1);
  make_vector(r->Ax,size);
  make_vector(r->Ai,size);
  return r;

error_malloc:
	freeSparseMatrix(r);
	return NULL;
}

SparseMatrix sparsem_pack(DOUBLE **m,unsigned int rows,unsigned cols)
{
  SparseMatrix r;
  if ((r = (SparseMatrix)malloc(sizeof(sparse_matrix))) == NULL)
    return NULL;
  r->rows = rows;
  r->cols = cols;

  unsigned int nz, count, i, j;
  nz = 0;
  for (i = 0;i < rows;i++)
    for (j = 0;j < cols;j++)
      if(m[i][j] != 0.0)
				nz++;
  
  r->nz = nz;
  r->size_of_Ax = nz;
  make_vector(r->Ax,nz);
  make_vector(r->Ap,cols+1);
  r->Ap[0] = 0;
  make_vector(r->Ai,nz);
	
  nz = 0;
  count = 0;
  
  for (j = 0;j < cols;j++)
	{
		for (i = 0;i < rows;i++)
			if (m[i][j] != 0.0)
			{
				r->Ax[count] = m[i][j];
				r->Ai[count] = i;
				count++;
				nz++;
			}
		r->Ap[j+1] = nz;
	}
  return r;

error_malloc:
	freeSparseMatrix(r);
	return NULL;
}

DOUBLE matrix_element(SparseMatrix r,unsigned int row,unsigned int col)
{
  unsigned int i;
  for (i = r->Ap[col]; i < r->Ap[col+1];i++)
    if(row == r->Ai[i])
      return r->Ax[i];
  return 0.0;
}

void print_sparsem_matrix(const char *fmt,SparseMatrix r)
{
  unsigned i, j;
  for (i = 0;i < r->rows;i++)
	{
		for (j = 0;j < r->cols;j++)
			printf(fmt,matrix_element(r,i,j));
		printf("\n");
	}
  printf("\n");  
}

void print_sparse_matrix_and_vector_to_txt_file(SparseMatrix r,DOUBLE *v,unsigned int length,const char *filename)
{
  unsigned i, j;
  FILE *fp;
  DOUBLE x;

  if ((fp = fopen(filename, "w")) == NULL)
  {
    fprintf(stderr, "cannot open file %s for writing\n",filename);
    return;
  }

  fprintf(fp,"A = {\n");
  for (i = 0;i < r->rows;i++)
  {
    fprintf(fp,"{");
    for (j = 0;j < r->cols;j++)
    {
      x = matrix_element(r,i,j);
      if (fabs(x) < 1.0E-5)
        x = 0.0;
      fprintf(fp,"%.12g",x);
      if (j < r->cols - 1)
        fprintf(fp,",");
    }
    fprintf(fp,"}");
    if (i < r->rows - 1)
        fprintf(fp,",");
    fprintf(fp,"\n");
  }
  fprintf(fp,"};\n\n");

  fprintf(fp,"B = {\n");
  for (i = 0;i < length;i++)
  {
    x = v[i];
    if (fabs(x) < 1.0E-4)
        x = 0.0;
    fprintf(fp,"%.12g",x);
    if (i < length - 1)
      fprintf(fp,",");
  }
  fprintf(fp,"};\n\n");
  fclose(fp);
}

void print_sparsem_matrix_elements(const char *fmt,SparseMatrix r)
{
  unsigned int i;
  printf("Rows: %u\n",r->rows);
  printf("Columns: %u\n",r->cols);
  printf("Number of non zero elements: %u\n",r->nz);
  if (fmt == NULL)
    return;
  
  printf("Ax: ");
  for (i = 0;i < r->nz;i++)
    printf(fmt,r->Ax[i]);
  printf("\nAp: ");
  for (i = 0;i <= r->cols;i++)
    printf("%u  ",r->Ap[i]);
  printf("\nAi: ");
  for (i = 0;i < r->nz;i++)
    printf("%u  ",r->Ai[i]);
  printf("\n");
}

static void my_write(void *ptr,size_t size,size_t nmemb,FILE *stream)
{
  if(fwrite(ptr,size,nmemb,stream) != nmemb)
  {
    fprintf(stderr,"Error writing a file\n");
    exit(EXIT_FAILURE);
  }
}

static void my_read(void *ptr,size_t size,size_t nmemb,FILE *stream)
{
  if(fread(ptr,size,nmemb,stream) != nmemb)
  {
    fprintf(stderr,"Error reading a file\n");
    exit(EXIT_FAILURE);
  }
}

void sparsem_write(FILE *fp,SparseMatrix r)
{
  my_write(&(r->rows),sizeof(unsigned int),1,fp);
  my_write(&(r->cols),sizeof(unsigned int),1,fp);
  my_write(&(r->nz),sizeof(unsigned int),1,fp);
	
  my_write(r->Ax,sizeof(DOUBLE),r->nz,fp);
  my_write(r->Ap,sizeof(unsigned int),r->cols+1,fp);
  my_write(r->Ai,sizeof(unsigned int),r->nz,fp);
}

SparseMatrix sparsem_read(FILE *fp)
{
  SparseMatrix r;
  if ((r = (SparseMatrix)malloc(sizeof(sparse_matrix))) == NULL)
    return NULL;
  my_read(&(r->rows),sizeof(unsigned int),1,fp);
  my_read(&(r->cols),sizeof(unsigned int),1,fp);
  my_read(&(r->nz),sizeof(unsigned int),1,fp);
  
  make_vector(r->Ax,r->nz);
  make_vector(r->Ap,r->cols+1);
  make_vector(r->Ai,r->nz);
  
  my_read(r->Ax,sizeof(DOUBLE),r->nz,fp);
  my_read(r->Ap,sizeof(unsigned int),r->cols+1,fp);
  my_read(r->Ai,sizeof(unsigned int),r->nz,fp);
  
  return r;
	
error_malloc:
	freeSparseMatrix(r);
	return NULL;
}

static int sparsem_append_element(SparseMatrix r,DOUBLE value,unsigned int row,unsigned int col)
/*
  The functions that uses sparsem_append_element  must append the elements in the compressed 
  sparse column order 
*/
{
  if(r->nz == r->size_of_Ax)
  {  
    expand_vector(r->Ax,r->nz + SIZEINCREMENT);
    expand_vector(r->Ai,r->nz + SIZEINCREMENT);
    r->size_of_Ax += SIZEINCREMENT;
  }
  r->Ax[r->nz] = value;
  r->Ai[r->nz] = row;
  r->nz++;
	return 1;

error_malloc:
	return 0;
}

int sparsem_change_element(SparseMatrix r,DOUBLE value,unsigned int row,unsigned int col,int mode)
{
  unsigned int pos, i;
  
  pos = r->Ap[col];
  while(pos < r->Ap[col+1])
  {
    if(row == r->Ai[pos])
    {
			if (mode == REPLACE)
				r->Ax[pos] = value;
			else
				r->Ax[pos] += value;
      return 1;
    }
    else if(row > r->Ai[pos])
      pos++;
    else
      break;
  }
  if(r->nz == r->size_of_Ax)
  {
    expand_vector(r->Ax,r->nz + SIZEINCREMENT);
    expand_vector(r->Ai,r->nz + SIZEINCREMENT);
    r->size_of_Ax += SIZEINCREMENT;
  }
  if(r->nz > pos)
  {
    memmove(r->Ai + pos + 1,r->Ai + pos,(r->nz - pos)*sizeof(unsigned int));
    memmove(r->Ax + pos + 1,r->Ax + pos,(r->nz - pos)*sizeof(DOUBLE));
  }
  r->Ax[pos] = value;
  r->Ai[pos] = row;
  r->nz++;

  for (i = col+1;i <= r->cols;i++)
    r->Ap[i]++;
	return 1;

error_malloc:
	return 0;
}

static void sparsem_add_to_element(SparseMatrix r,DOUBLE value,unsigned int row,unsigned int col)
{
	unsigned int pos;
 
	pos = r->Ap[col];
	while(pos < r->Ap[col+1])
	{
		if(row == r->Ai[pos])
		{
			r->Ax[pos] += value;
			return;
		}
		pos++;
	}
}

DOUBLE *sparsem_multiply_by_vector(SparseMatrix r,DOUBLE *x)
{
  DOUBLE *result;
  unsigned int i, j;
  
  make_vector(result,r->cols);
  for (i = 0;i < r->cols;i++)
    for (j = r->Ap[i];j < r->Ap[i+1];j++)
      result[r->Ai[j]] += r->Ax[j]*x[i];
  return result;
	
error_malloc:
	return NULL;
}

void free_sparsem_matrix(SparseMatrix *r)
{
  if(*r == NULL)
    return;
  free_vector((*r)->Ap);
  free_vector((*r)->Ai);
  free_vector((*r)->Ax);
  free(*r);
  *r = NULL;
}

TripletForm triplet_form_init(unsigned int size)
{
  TripletForm t;
  if ((t = (TripletForm)malloc(sizeof(triplet_form))) == NULL)
    return NULL;
  make_vector(t->T,size);
  t->size_of_T = size;
  t->elements = 0;
  t->rows = 0;
  t->cols = 0;
  return t;

error_malloc:
	freeTripletForm(t);
	return NULL;
}

int triplet_form_append_element(TripletForm t,DOUBLE value,unsigned int row,unsigned int col)
{
  if(t->elements == t->size_of_T)
  {  
    expand_vector(t->T,t->elements + SIZEINCREMENT);
    t->size_of_T += SIZEINCREMENT;
  }
  t->T[t->elements].row = row;
  t->T[t->elements].col = col;
  t->T[t->elements].value = value;
  t->elements++;

  if(++row > t->rows)
    t->rows = row;
  if(++col > t->cols)
    t->cols = col;
	return 1;

error_malloc:
	return 0;
}


TripletForm triplet_form_read(const char *filename,unsigned char symmetric,unsigned int zero_based)
{
  FILE *fp;
  TripletForm t;
  unsigned int i, j, m;
  DOUBLE x;
  
  t = NULL;
  fp = NULL;
  m = (zero_based == 0) ? 1 : 0; 

  if ((fp = fopen(filename,"r")) == NULL)
    goto final;
  if ((t = triplet_form_init(SIZEINCREMENT)) == NULL)
    goto final;

  while (fscanf(fp,"%u %u %lg\n",&i,&j,&x) == 3)
	{
		if (! triplet_form_append_element(t,x,i - m,j - m))
			goto error;
		if ((i != j) && (symmetric != 0))
			if (! triplet_form_append_element(t,x,j - m,i - m))
				goto error;
	}
  
final:
  if (fp != NULL)
    fclose(fp);
  return t;

error:
	freeTripletForm(t);
	return NULL;
}

static int smcmpfunc(const void *a,const void *b)
{
  MatrixElement *da = (MatrixElement *)a;
  MatrixElement *db = (MatrixElement *)b;
  if(da->col > db->col)
    return 1;
  if(da->col < db->col)
    return -1;
  if(da->row > db->row)
    return 1;
  if(da->row < db->row)
    return -1;
  return 0;
}

SparseMatrix sparsem_from_triplet_form(TripletForm t,int mode)
{
  SparseMatrix r;
  unsigned int col, row, n, first;
  qsort(t->T,t->elements,sizeof(MatrixElement),smcmpfunc);
  if ((r = sparsem_init(t->rows,t->cols,t->size_of_T)) == NULL)
    return NULL;
  n = 0;
  for (col = 0;col < t->cols;col++) 
	{
		while(t->T[n].col == col)
		{
			DOUBLE x = 0.0;
			row = t->T[n].row;
			while ((t->T[n].col == col) && (t->T[n].row == row))
			{
				if (mode == REPLACE)
					x = t->T[n].value;
				else
					x += t->T[n].value;
				n++;
			}
			if (! sparsem_append_element(r,x,row,col))
					goto error;
		}
		r->Ap[col+1] = r->nz;
	}
  return r;

error:
	freeSparseMatrix(r);
	return NULL;
}

void print_triplet_form(const char *fmt,TripletForm t)
{
  unsigned i;
  printf("Elements: %d    Rows: %d      Columns: %d\n",t->elements,t->rows,t->cols);
  for (i=0;i < t->elements;i++)
	{
		printf("Row: %d    Column: %d      ",t->T[i].row,t->T[i].col);
		printf(fmt,t->T[i].value);
	}
  printf("\n");
}

void free_triplet_form(TripletForm *t)
{
  if(*t == NULL)
    return;
  free_vector((*t)->T);
  free(*t);
  *t = NULL;
}

DOUBLE *umfpack_solve(SparseMatrix r,DOUBLE *b)
{
  void *Symbolic, *Numeric;
  DOUBLE Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];
  int status;
  DOUBLE *x;

  x = NULL;
  status = umfpack_di_symbolic (r->rows,r->cols,r->Ap,r->Ai,r->Ax,&Symbolic,Control,Info);
  if (status != UMFPACK_OK)
    goto error_malloc;
  
  status = umfpack_di_numeric (r->Ap,r->Ai,r->Ax,Symbolic,&Numeric,Control,Info);
  umfpack_di_free_symbolic (&Symbolic);
  if (status != UMFPACK_OK)
		goto error_malloc;
  
  make_vector(x,r->cols);
  status = umfpack_di_solve (UMFPACK_A,r->Ap,r->Ai,r->Ax,x,b,Numeric,Control,Info) ;
  umfpack_di_free_numeric(&Numeric);
  if (status != UMFPACK_OK)
		goto error_malloc;

  return x;

error_malloc:
  printf("UMFPack error code: %d\n",status);
	free_vector(x);
	return NULL;
}

