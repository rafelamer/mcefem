
/**************************************************************************************
* Filename:   lagrange.c
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
#include <stdlib.h>

/*
	One dimension Lagrange polynomials in the interval [0,1]
*/

DOUBLE Lagrange1D(DOUBLE x,unsigned int i,unsigned int n)
	/*
		i, k are unsigned int and we have take care with the subtraction i - k
		If we put 1.0 * i - 1.0 * k, they are converted to floats.
	*/
{
	unsigned int k;
	DOUBLE r = 1.0;
	
	if (i > n)
		return NAN;
  
	for (k = 0;k <= n;k++)
	{
		if (k == i)
			continue;
		r *= (n * x - k) / (1.0 * i - 1.0 * k);
	}
	return r;
}

DOUBLE DerivativeLagrange1D(DOUBLE x,unsigned int i,unsigned int n)
{
	unsigned int k, j;
	DOUBLE s, r;

	if (i > n)
		return NAN;
  
	s = 0.0;
	for (k = 0;k <= n;k++)
	{
		if(k == i)
			continue;
		r = 1.0;
		for (j = 0;j <= n;j++)
		{
			if ((j == k) || (j == i))
				continue;
			r *= (n * x - j) / (1.0 * i - 1.0 * j);
		}
		s += r / (1.0 * i - 1.0 * k);
	}
  return n * s;
}

Lagrange1DValues LagrangeAtGaussPoints(int degree,GaussQuadrature qdat,int npoints)
{
	Lagrange1DValues lg;
	unsigned int k, i;

	lg = NULL;
	make_matrix(lg,degree + 1,npoints);

	for (k = 0;k <= degree;k++)
	{
		for (i = 0;i < npoints;i++)
		{
			DOUBLE x = (1 + qdat[i].x) / 2;
			lg[k][i].value =  Lagrange1D(x,k,degree);
			lg[k][i].dvalue =  DerivativeLagrange1D(x,k,degree);
		}
	}
	return lg;
	
error_malloc:
	free_matrix(lg);
	return NULL;
}

Lagrange1DValues LagrangeAtPointsInBasicInterval(int degree,int npoints)
{
	Lagrange1DValues lg;
	unsigned int i, j;
	DOUBLE h;

	lg = NULL;
	make_matrix(lg,degree + 1,npoints);
	h = 1.0 / (npoints - 1);

	for (i = 0;i <= degree;i++)
	{
		for (j = 0;j < npoints;j++)
		{
			DOUBLE x = j * h;		
			lg[i][j].value =  Lagrange1D(x,i,degree);
			lg[i][j].dvalue =  DerivativeLagrange1D(x,i,degree);
		}
	}
	return lg;
	
error_malloc:
	free_matrix(lg);
	return NULL;
}

/*
	Two dimension Lagrange polynomials in the basic triangle
*/

#define F(t)  for (k = 0;k < i;k++)				                            \
	t *= (n * x - k) / (1.0 * i - 1.0 * k);

#define G(t)  for (k = 0;k < j;k++)				                            \
	t *= (n * y - k) / (1.0 * j - 1.0 * k);

#define H(t)  for (k = i + j + 1;k <= n;k++)                                \
	t *= (n * x + n * y - k) / (1.0 * i + 1.0 * j - 1.0 * k);

#define DFGXY(t,m)    s = 0.0;                                              \
	for (k = 0;k < m;k++)				                                    \
	{                                                                       \
		r = 1.0;                                        					\
		for (p = 0;p < m;p++)                                               \
		{                                                                   \
			if (p == k)                                                     \
				continue;                                                   \
			r *= (n * t - p) / (1.0 * m - 1.0 * p);                         \
		}                                                                   \
		s += r / (1.0 * m - 1.0 * k);                                       \
	}

#define DHXY    s = 0.0;                                                    \
	for (k = i + j + 1;k <= n;k++)			                                \
	{                                                                       \
		r = 1.0;                                                            \
		for (p = i + j + 1;p <= n;p++)                                      \
		{                                                                   \
			if (p == k)                                                     \
				continue;                                                   \
			r *= (n * x + n * y  - p) / (1.0 * i + 1.0 * j - 1.0 * p);      \
		}                                                                   \
		s += r / (1.0 * i + 1.0 * j - 1.0 * k);                             \
	}

DOUBLE Lagrange2D(DOUBLE x,DOUBLE y,unsigned int i,unsigned j,unsigned int n)
{
	unsigned int k;
	DOUBLE r = 1.0;

	if ((i > n) || (j > n - i))
    	return NAN;

	F(r);
	G(r);
	H(r);
	return r;
}

DOUBLE DerivativeXLagrange2D(DOUBLE x,DOUBLE y,unsigned int i,unsigned j,unsigned int n)
{
	DOUBLE r, s, result;
	unsigned int k, p;
  
	if ((i > n) || (j > n - i))
		return NAN;

	DFGXY(x,i);
	G(s);
	H(s);
	result = s;

	DHXY;
	F(s);
	G(s);

	result += s;
	return n * result;
}

DOUBLE DerivativeYLagrange2D(DOUBLE x,DOUBLE y,unsigned int i,unsigned j,unsigned int n)
{
	DOUBLE r, s, result;
	unsigned int k, p;
  
	if ((i > n) || (j > n - i))
		return NAN;

	DFGXY(y,j);
	F(s);
	H(s);
	result = s;
  
	DHXY;
	F(s);
	G(s);

	result += s;
	return n * result;
}

Lagrange2DValues LagrangeAtTWBPoints(int elements,TWBGaussQuadrature qdat,int npoints)
{
	int i, j, k;
	Lagrange2DValues lv;

	make_triangular_matrix(lv,elements + 1);
	for (i = 0; i <= elements;i++)
	{
		for (j = 0;j <= elements - i;j++)
		{
			make_vector(lv[i][j],npoints);
			for (k = 0;k < npoints;k++)
			{
				DOUBLE x, y;
				x = qdat[k].lambda1;
				y = qdat[k].lambda2;
				lv[i][j][k].value =  Lagrange2D(x,y,i,j,elements); 
				lv[i][j][k].dxvalue = DerivativeXLagrange2D(x,y,i,j,elements); 
				lv[i][j][k].dyvalue = DerivativeYLagrange2D(x,y,i,j,elements); 
			}
		}
	}
	return lv;

error_malloc:
	freeLagrange2DValues(lv,elements);
	return NULL;
}

Lagrange2DValues LagrangeAtPointsInBasicTriangle(int degree,int npoints)
{
	int i, j, k, l, values;
	Lagrange2DValues lv;
	DOUBLE h;

	values = npoints * (npoints + 1) / 2;
	make_triangular_matrix(lv,degree + 1);
	h = 1.0 / (npoints - 1);
	for (i = 0; i <= degree;i++)
	{
		for (j = 0;j <= degree - i;j++)
		{
			make_vector(lv[i][j],values);
			int count = 0;
			for (k = 0;k < npoints;k++)
			{
				for (l = 0;l < npoints - k;l++)
				{
					DOUBLE x, y;
					x = k * h;
					y = l * h;
					lv[i][j][count].value =  Lagrange2D(x,y,i,j,degree);
					lv[i][j][count].dxvalue = DerivativeXLagrange2D(x,y,i,j,degree); 
					lv[i][j][count].dyvalue = DerivativeYLagrange2D(x,y,i,j,degree);
					count++;
				} 
			}
		}
	}
	return lv;

error_malloc:
	freeLagrange2DValues(lv,degree);
	return NULL;
}

void free_lagrange_1d_values(Lagrange1DValues *data)
{
	free_matrix(*data);
}

void free_lagrange_2d_values(Lagrange2DValues *data,int elements)
{
	int i, j;
	for (i = 0; i <= elements;i++)
		for (j = 0;j <= elements - i;j++)
			free_vector((*data)[i][j]);
	free_matrix(*data);
}

