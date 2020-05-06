/**************************************************************************************
* Filename:   tfgfem.h
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
*             See https://www.gnu.org/licenses/
***************************************************************************************/
#include <tfgfem.h>

DOUBLE gauss_integrate_function(DOUBLE a,DOUBLE b,FunctionX f,int points)
{
	GaussQuadrature qdat;
	DOUBLE sum, alfa, beta;

	qdat = gauss_data(&points);
	sum = 0.0;
	alfa = 0.5 * (b - a);
	beta = 0.5 * (b + a);
	while(qdat->w != -1)
	{
		sum += qdat->w * f(alfa * qdat->x + beta);
		qdat++;
	}
	return alfa * sum;
}

DOUBLE gauss_integrate_evaluator(DOUBLE a,DOUBLE b,void *evaluator,int points)
{
	GaussQuadrature qdat;
	DOUBLE sum, alfa, beta;

	qdat = gauss_data(&points);
	sum = 0.0;
	alfa = 0.5 * (b - a);
	beta = 0.5 * (b + a);
	while(qdat->w != -1)
	{
		sum += qdat->w * evaluator_evaluate_x(evaluator,alfa * qdat->x + beta);
		qdat++;
	}
	return alfa * sum;
}

DOUBLE twb_integrate_triangle_function(PTriangle tr,FunctionXY f,TWBGaussQuadrature qdat)
{
 	DOUBLE x[3], y[3];
	DOUBLE l1, l2, X, Y, sum;
	int i;

	sum = 0.0;
	for (i = 0; i < 3; i++) 
	{
		x[i] = tr->n[i]->x;
		y[i] = tr->n[i]->y;
	}
  
	while (qdat->weight != -1) 
	{
		l1 = qdat->lambda1;
		l2 = qdat->lambda2;
		X = x[0] + (x[1] - x[0]) * l1 + (x[2] - x[0]) * l2;
		Y = y[0] + (y[1] - y[0]) * l1 + (y[2] - y[0]) * l2;
		sum += qdat->weight * f(X,Y);
		qdat++;
	}
	return 0.5 * tr->area * sum;
}

DOUBLE twb_integrate_triangle_evaluator(PTriangle tr,void *evaluator,TWBGaussQuadrature qdat)
{
	DOUBLE x[3], y[3];
	DOUBLE l1, l2, X, Y, sum;
	int i;

	sum = 0.0;
	for (i = 0; i < 3; i++) 
	{
		x[i] = tr->n[i]->x;
		y[i] = tr->n[i]->y;
	}
  
	while (qdat->weight != -1) 
	{
		l1 = qdat->lambda1;
		l2 = qdat->lambda2;
		X = x[0] + (x[1] - x[0]) * l1 + (x[2] - x[0]) * l2;
		Y = y[0] + (y[1] - y[0]) * l1 + (y[2] - y[0]) * l2;
		sum += qdat->weight * evaluator_evaluate_x_y(evaluator,X,Y);
		qdat++;
	}

	return 0.5 * tr->area * sum;
}


DOUBLE twb_integrate_mesh_function(DataMesh mesh,FunctionXY f,int degree)
{
	DOUBLE sum = 0.0;
	int i, d;
	TWBGaussQuadrature qdat;

	qdat = get_twb_data(&degree,NULL);
	for (i = 0;i < mesh->ntriangles;i++)
		sum += twb_integrate_triangle_function(&mesh->triangles[i],f,qdat);
	return sum;
}

DOUBLE twb_integrate_mesh_evaluator(DataMesh mesh,void *evaluator,int degree)
{
	DOUBLE sum = 0.0;
	int i, d;
	TWBGaussQuadrature qdat;

	qdat = get_twb_data(&degree,NULL);
	for (i = 0;i < mesh->ntriangles;i++)
		sum += twb_integrate_triangle_evaluator(&mesh->triangles[i],evaluator,qdat);
	return sum;
}

#define Value(m,n) (lv2[(m)][(n)][count].value)
#define DXValue(m,n) (a22 * lv2[(m)][(n)][count].dxvalue - a21 * lv2[(m)][(n)][count].dyvalue)
#define DYValue(m,n) (-a12 * lv2[(m)][(n)][count].dxvalue + a11 * lv2[(m)][(n)][count].dyvalue)

#define LagrangeOneTerm(m,n,derivative)    \
	switch ((derivative))                  \
	{                                      \
		case NODERIVATIVE:                 \
		{                                  \
			t *= Value(m,n);               \
			break;                         \
		}                                  \
		case DXDERIVATIVE:                 \
		{                                  \
			t *= DXValue(m,n);             \
			break;                         \
		}                                  \
		case DYDERIVATIVE:                 \
		{                                  \
			t *= DYValue(m,n);             \
			break;                         \
		}                                  \
		default:                           \
			t = NAN;                       \
	}

#define LagrangeTwoTerms(m1,n1,m2,n2,derivatives)      \
	switch ((derivatives))                             \
	{                                                  \
		case NODERIVATIVES:                            \
		{                                              \
			t *= Value(m1,n1) * Value(m2,n2);          \
			break;                                     \
		}                                              \
		case DXNODERIVATIVES:                          \
		{                                              \
			t *= DXValue(m1,n1) * Value(m2,n2);        \
			break;                                     \
		}                                              \
		case NODXDERIVATIVES:                          \
		{                                              \
			t *= Value(m1,n1) * DXValue(m2,n2);        \
			break;                                     \
		}                                              \
		case DXDXDERIVATIVES:                          \
		{                                              \
			t *= DXValue(m1,n1) * DXValue(m2,n2);      \
			break;                                     \
		}                                              \
		case DYNODERIVATIVES:                          \
		{                                              \
			t *= DYValue(m1,n1) * Value(m2,n2);        \
			break;                                     \
		}                                              \
		case NODYDERIVATIVES:                          \
		{                                              \
			t *= Value(m1,n1) * DYValue(m2,n2);        \
			break;                                     \
		}                                              \
		case DYDYDERIVATIVES:                          \
		{                                              \
			t *= DYValue(m1,n1) * DYValue(m2,n2);      \
			break;                                     \
		}                                              \
		case DXDYDERIVATIVES:                          \
		{                                              \
			t *= DXValue(m1,n1) * DYValue(m2,n2);      \
			break;                                     \
		}                                              \
		case DYDXDERIVATIVES:                          \
		{                                              \
			t *= DYValue(m1,n1) * DXValue(m2,n2);      \
			break;                                     \
		}                                              \
		case NABLADERIVATIVES:                         \
		{                                              \
			t *= (DXValue(m1,n1) * DXValue(m2,n2)      \
				+ DYValue(m1,n1) * DYValue(m2,n2));    \
			break;                                     \
		}                                              \
		default:                                       \
			t = NAN;                                   \
	}

#define integral_over_triangle_one_lagrange(fun,m1,n1,d1,res)  {                       \
	DOUBLE x[3], y[3];                                                                 \
	DOUBLE la1, la2, X, Y, sum, a11, a12, a21, a22, det, factor;                       \
	int z, count;                                                                      \
	sum = 0.0;                                                                         \
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
		t = q2->weight * runFunctionTable2D(fun,X,Y);                                  \
		LagrangeOneTerm(m1,n1,d1);                                                     \
		sum += t;                                                                      \
		q2++;                                                                          \
		count++;                                                                       \
	}                                                                                  \
	if (d1 == 0)                                                                       \
		factor = tr->area;                                                             \
	else                                                                               \
		factor = 0.5 * ((det > 0) - (det < 0));                                        \
	res = 0.5 * factor * sum;                                                          \
}

#define integral_over_triangle_two_lagrange(fun1,fun2,m1,n1,m2,n2,derivative,res)  {   \
	DOUBLE x[3], y[3];                                                                 \
	DOUBLE la1, la2, X, Y, sum, a11, a12, a21, a22, det, factor;                       \
	int z, count;                                                                      \
	sum = 0.0;                                                                         \
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
		t = q2->weight * runFunctionTable2D(fun1,X,Y);                                 \
		LagrangeTwoTerms(m1,n1,m2,n2,derivative);                                      \
		sum += t;                                                                      \
        if (fun2 != NULL)                                                              \
        {                                                                              \
        	t = q2->weight * runFunctionTable2D(fun2,X,Y);                             \
        	LagrangeTwoTerms(m1,n1,m2,n2,NODYDERIVATIVES);                             \
			sum += t;                                                                  \
        }                                                                              \
		q2++;                                                                          \
		count++;                                                                       \
	}                                                                                  \
	if (derivative == NODERIVATIVES)                                                   \
		factor = tr->area;                                                             \
	else if (derivative <= NODYDERIVATIVES)                                            \
		factor = 0.5 * ((det > 0) - (det < 0));                                        \
	else                                                                               \
		factor = 0.25/tr->area;                                                        \
	res = 0.5 * factor * sum;                                                          \
}

DOUBLE integralOverOneTriangleWithOneLagrange(PTriangle tr,unsigned int m1,unsigned int n1,unsigned int d1,
											  FunctionTable item,TWBGaussQuadrature q2,Lagrange2DValues lv2)
{
	DOUBLE res;
	integral_over_triangle_one_lagrange(item,m1,n1,d1,res);
	return res;
}

DOUBLE integralOverOneTriangleWithTwoLagrange(PTriangle tr,unsigned int m1,unsigned int n1,
											  unsigned int m2,unsigned int n2,unsigned int derivative,  
											  FunctionTable item,TWBGaussQuadrature q2,Lagrange2DValues lv2)
{
	DOUBLE res;
	integral_over_triangle_two_lagrange(item,NULL,m1,n1,m2,n2,derivative,res);
	return res;
}

DOUBLE integralOverOneTriangleWithTwoLagrangeTest(PTriangle tr,unsigned int m1,unsigned int n1,
												  unsigned int m2,unsigned int n2,
												  FunctionTable item,TWBGaussQuadrature q2,Lagrange2DValues lv2)
{
	DOUBLE res;
	FunctionTable next;
	next = item->hh.next;
	integral_over_triangle_two_lagrange(item,next,m1,n1,m2,n2,NODXDERIVATIVES,res);
	return res;
}

#define integral_function_over_one_edge(item,edge,r)  {                   \
	DOUBLE X, Y;                                                          \
	r = 0.0;                                                              \
	while(q1->w != -1)                                                    \
	{                                                                     \
		DOUBLE t;                                                         \
		t = 0.5 * q1->x + 0.5;                                            \
		X = tr->e[edge]->n[0]->x + t * tr->ex[edge];                      \
		Y = tr->e[edge]->n[0]->y + t * tr->ey[edge];                      \
		r += q1->w * runFunctionTable2D(item,X,Y);                        \
		q1++;                                                             \
	}                                                                     \
	r *= 0.5 *sqrt(tr->ex[edge]*tr->ex[edge]+tr->ey[edge]*tr->ey[edge]);  \
}

DOUBLE integralFunctionOverOneEdge(PTriangle tr,int edge,FunctionTable item,GaussQuadrature q1)
{
	DOUBLE r;
	integral_function_over_one_edge(item,edge,r);
	printf("First node of the edge (%.16g ,%.16g)\n",tr->e[edge]->n[0]->x,tr->e[edge]->n[0]->y);
	printf("Second node of the edge (%.16g ,%.16g)\n",tr->e[edge]->n[1]->x,tr->e[edge]->n[1]->y);
	printf("Computed second node of the edge (%.16g ,%.16g)\n",tr->e[edge]->n[0]->x + tr->ex[edge],tr->e[edge]->n[0]->y + tr->ey[edge]);
	return r;
}

#define integral_over_one_edge_two_lagrange(item,edge,m1,n1,r) {                       \
	DOUBLE X, Y;                                                                       \
	int n = 0;                                                                         \
	r = 0.0;                                                                           \
	while(q1->w != -1)                                                                 \
	{                                                                                  \
		DOUBLE t;                                                                      \
		t = 0.5 * q1->x + 0.5;                                                         \
		X = tr->e[edge]->n[0]->x + t * tr->ex[edge];                                   \
		Y = tr->e[edge]->n[0]->y + t * tr->ey[edge];                                   \
		r += q1->w * runFunctionTable2D(item,X,Y) * lv[m1][n].value * lv[n1][n].value; \
		q1++;                                                                          \
		n++;                                                                           \
	}                                                                                  \
	r *= 0.5 *sqrt(tr->ex[edge]*tr->ex[edge]+tr->ey[edge]*tr->ey[edge]);               \
}

DOUBLE integralOverOneEdgeTwoLagrange(PTriangle tr,int edge,unsigned int m1,unsigned int n1,
									  FunctionTable item,GaussQuadrature q1,Lagrange1DValues lv)
{
	DOUBLE r;
	integral_over_one_edge_two_lagrange(item,edge,m1,n1,r);
	printf("First node of the edge (%.16g ,%.16g)\n",tr->e[edge]->n[0]->x,tr->e[edge]->n[0]->y);
	printf("Second node of the edge (%.16g ,%.16g)\n",tr->e[edge]->n[1]->x,tr->e[edge]->n[1]->y);
	printf("Computed second node of the edge (%.16g ,%.16g)\n",tr->e[edge]->n[0]->x + tr->ex[edge],tr->e[edge]->n[0]->y + tr->ey[edge]);
	return r;
}
