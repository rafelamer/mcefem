/**************************************************************************************
* Filename:   regions.c
* Author:     Rafel Amer (rafel.amer AT upc.edu)
* Copyright:  Rafel Amer 2018
* Disclaimer: This code is presented "as is" and it has been written to 
*             implement the Finite Element Method in dimension 2. It
*             has been writen educational purposes.
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
*	      See https://www.gnu.org/licenses/
***************************************************************************************/
#include <tfgfem.h>
#include <math.h>
#include <time.h>
#include <matheval.h>

#define PI 4*atan(1)
#define AREEQUAL(a,b) !xmlStrcmp((a),(const xmlChar *)(b))

DataRegion data_region_new()
{
	DataRegion r;
	if ((r = (DataRegion)malloc(sizeof(data_region))) == NULL)
		return NULL;
	return r;
}

void free_data_region(DataRegion *r)
{
	if (*r == NULL)
		return;
	free_vector((*r)->points);
	free_vector((*r)->segments);
	free_vector((*r)->holes);
	free(*r);
	*r = NULL;
}

static void parsePoints(xmlNodePtr node,DataRegion r,int bc)
{
	xmlNodePtr cur;
	int i, n, m, found, alloc, size;
	data_point *points;
	data_segment *edges;

	n = r->npoints;
	m = r->nsegments;
	alloc = 256;
	make_vector(points,alloc);
	i = 0;
	cur = node->xmlChildrenNode;
	while (cur != NULL)
	{
		if (AREEQUAL(cur->name,"point"))
		{
			DOUBLE x, y;
			int pbc, sbc;
			x = get_double_value(cur,"x",&found,XML_VALUE_REQUIRED);    
			y = get_double_value(cur,"y",&found,XML_VALUE_REQUIRED);
			pbc = get_int_value(cur,"bc",&found,XML_VALUE_OPTIONAL);
			if (found == 0)
				pbc = bc;
			sbc = get_int_value(cur,"sbc",&found,XML_VALUE_OPTIONAL);
			if (found == 0)
				sbc = bc;
			if (i == alloc)
	    	{
				alloc += 256;
				expand_vector(points,alloc);
	    	}
			points[i].point_no = n + i;
			points[i].x = x;
			points[i].y = y;
			points[i].bc = pbc;
			points[i].sbc = sbc;
			i++;
		}
		cur = cur->next;
	}
	
	/*
    	Append edges
	*/
	size = i;
	r->npoints = n + size;  
	make_vector(edges,size);
	
	for (i = 0;i < size;i++)
	{
		edges[i].segment_no = m + i;
		edges[i].point_no_1 = n + i;
		edges[i].point_no_2 = n + ((i + 1) % size);
		edges[i].bc = points[i].sbc;
	}
	append_to_vector(r->points,points,n,size,FREE_APPENDED);
	append_to_vector(r->segments,edges,m,size,FREE_APPENDED);
	r->nsegments = m + size;
	return;

error_malloc:
	fprintf(stderr,"Error in a malloc or calloc call\n");
	exit(EXIT_FAILURE);
}

static void parsePolygon(xmlNodePtr node,DataRegion r)
{
	xmlChar *hole;
	int bc, found;
  
	hole =  get_string_value(node,"ishole",&found,XML_VALUE_OPTIONAL);
	if (found && AREEQUAL(hole,"true"))
	{
		data_hole *h;
		double x, y;
		int n;
		n = r->nholes;
		x = get_double_value(node,"x",&found,XML_VALUE_REQUIRED);
		y = get_double_value(node,"y",&found,XML_VALUE_REQUIRED);
		make_vector(h,1);
		h->x = x;
		h->y = y;
		append_to_vector(r->holes,h,n,1,FREE_APPENDED);
		r->nholes = n+1;
	}

	bc = get_int_value(node,"bc",&found,XML_VALUE_OPTIONAL);
	if (found == 0)
		bc = FEM_BC_DIRICHLET;
  
	parsePoints(node,r,bc);
	return;

error_malloc:
	fprintf(stderr,"Error in a malloc or calloc call\n");
	exit(EXIT_FAILURE);	
}

static void parseCircle(xmlNodePtr node,DataRegion r)
{
	xmlChar *hole;
	xmlNodePtr cur;
	int found, number, bc;
  
	hole =  get_string_value(node,"ishole",&found,XML_VALUE_OPTIONAL);
	if (found && AREEQUAL(hole,"true"))
	{
		data_hole *h;
		double x, y;
		int n;
		n = r->nholes;
		x = get_double_value(node,"x",&found,XML_VALUE_REQUIRED);
		y = get_double_value(node,"y",&found,XML_VALUE_REQUIRED);
		make_vector(h,1);
		h->x = x;
		h->y = y;
		append_to_vector(r->holes,h,n,1,FREE_APPENDED);
		r->nholes = n+1;
	}
	bc = get_int_value(node,"bc",&found,XML_VALUE_OPTIONAL);
	if (found == 0)
		bc = FEM_BC_DIRICHLET;
	
	number = 0;
	cur =  node->xmlChildrenNode;
	while (cur != NULL)
	{
		if (AREEQUAL(cur->name,"center"))
		{
			int segments, n, m, i;
			double x, y, radius;
			data_segment *edges;
			data_point *points;
			
			number++;
			if (number > 1)
			{
				fprintf(stderr,"Only one center is allowed in regions of type circle\n");
				exit(EXIT_FAILURE); 
			}
			
			n = r->npoints;
			m = r->nsegments;
			x = get_double_value(cur,"x",&found,XML_VALUE_REQUIRED);     
			y = get_double_value(cur,"y",&found,XML_VALUE_REQUIRED);
			radius = get_double_value(cur,"radius",&found,XML_VALUE_REQUIRED);
			segments = get_int_value(cur,"segments",&found,XML_VALUE_REQUIRED);
			if(segments < 3)
				segments = 3;
			make_vector(edges,segments);
			make_vector(points,segments);
			for (i = 0;i < segments;i++)
			{
				DOUBLE t = 2*i*PI/segments;
				double p,q;
				points[i].point_no = n+i;
				p = radius*cos(t);
				q = radius*sin(t);	
				points[i].x = x + p;
				points[i].y = y + q;
				points[i].bc = bc;
				edges[i].segment_no = m+i;
				edges[i].point_no_1 = n+i;
				edges[i].point_no_2 = n+(i+1) % segments;
				edges[i].bc = bc;
			}
			append_to_vector(r->points,points,n,segments,FREE_APPENDED);
			append_to_vector(r->segments,edges,m,segments,FREE_APPENDED);
			r->npoints = segments + n;
			r->nsegments = segments + m;
		}
		cur = cur->next;
	}
	if (number == 0)
	{
		fprintf(stderr,"One center is required in regions of type circle\n");
		exit(EXIT_FAILURE); 
	}
	return;

error_malloc:
	fprintf(stderr,"Error in a malloc or calloc call\n");
	exit(EXIT_FAILURE);
}

static void parseSquare(xmlNodePtr node,DataRegion r)
{
	xmlChar *hole;
	xmlNodePtr cur;
	int found, number, bc;
  
	hole =  get_string_value(node,"ishole",&found,XML_VALUE_OPTIONAL);
	if (found && AREEQUAL(hole,"true"))
	{
		data_hole *h;
		double x, y;
		int n;
		n = r->nholes;
		x = get_double_value(node,"x",&found,XML_VALUE_REQUIRED);
		y = get_double_value(node,"y",&found,XML_VALUE_REQUIRED);
		make_vector(h,1);
		h->x = x;
		h->y = y;
		append_to_vector(r->holes,h,n,1,FREE_APPENDED);
		r->nholes = n+1;
	}
	bc = get_int_value(node,"bc",&found,XML_VALUE_OPTIONAL);
	if (found == 0)
	bc = FEM_BC_DIRICHLET;
	
	number = 0;
	cur =  node->xmlChildrenNode;
	while (cur != NULL)
	{
		if (AREEQUAL(cur->name,"center"))
		{
			int n, m, i;
			DOUBLE x, y, side, rotate, c, s, p, q;
			data_segment *edges;
			data_point *points;
			
			number++;
			if (number > 1)
		    {
				fprintf(stderr,"Only one center is allowed in regions of type square\n");
				exit(EXIT_FAILURE); 
			}
			n = r->npoints;
			m = r->nsegments;
			x = get_double_value(cur,"x",&found,XML_VALUE_REQUIRED);     
			y = get_double_value(cur,"y",&found,XML_VALUE_REQUIRED);
			side = get_double_value(cur,"side",&found,XML_VALUE_REQUIRED);
			rotate = get_double_value(cur,"rotate",&found,XML_VALUE_OPTIONAL);
			if (found == 0)
				rotate = 0.0;
			make_vector(edges,4);
			make_vector(points,4);
			c = cos(PI*rotate/180);
			s = sin(PI*rotate/180);
			/*
				First point
			 */
			points[0].point_no = n;
			p = side/2.0;
			q = -side/2.0;
			points[0].x = x+c*p-s*q;
			points[0].y = y+s*p+c*q;
			points[0].bc = bc;
			/*
				Second point
			 */
			points[1].point_no = n + 1;
			p = side/2.0;
			q = side/2.0;
			points[1].x = x+c*p-s*q;
			points[1].y = y+s*p+c*q;
			points[1].bc = bc;
			/*
				Third point
			 */
			points[2].point_no = n + 2;
			p = -side/2.0;
			q = side/2.0;
			points[2].x = x+c*p-s*q;
			points[2].y = y+s*p+c*q;
			points[2].bc = bc;
			/*
				Fourth point
			 */
			points[3].point_no = n + 3;
			p = -side/2.0;
			q = -side/2.0;
			points[3].x = x+c*p-s*q;
			points[3].y = y+s*p+c*q;
			points[3].bc = bc;
			for (i = 0;i < 4;i++)
			{
				edges[i].segment_no = m+i;
				edges[i].point_no_1 = n+i;
				edges[i].point_no_2 = n + (i+1) % 4;
				edges[i].bc = bc;
			}
			append_to_vector(r->points,points,n,4,FREE_APPENDED);
			append_to_vector(r->segments,edges,m,4,FREE_APPENDED);
			r->npoints = 4 + n;
			r->nsegments = 4 + m;
		}
		cur = cur->next;
	}	
	if (number == 0)
	{
		fprintf(stderr,"One center is required in regions of type square\n");
		exit(EXIT_FAILURE); 
	}
	return;

error_malloc:
	fprintf(stderr,"Error in a malloc or calloc call\n");
	exit(EXIT_FAILURE);
}

static void parseEllipse(xmlNodePtr node,DataRegion r)
{
	xmlChar *hole;
	xmlNodePtr cur;
	int found, number, bc;
  
	hole =  get_string_value(node,"ishole",&found,XML_VALUE_OPTIONAL);
	if (found && AREEQUAL(hole,"true"))
	{
		data_hole *h;
		DOUBLE x, y;
		int n;
		n = r->nholes;
		x = get_double_value(node,"x",&found,XML_VALUE_REQUIRED);
		y = get_double_value(node,"y",&found,XML_VALUE_REQUIRED);
		make_vector(h,1);
		h->x = x;
		h->y = y;
		append_to_vector(r->holes,h,n,1,FREE_APPENDED);
 		r->nholes = n+1;
	}
	bc = get_int_value(node,"bc",&found,XML_VALUE_OPTIONAL);
	if (found == 0)
		bc = FEM_BC_DIRICHLET;

	number = 0;
	cur =  node->xmlChildrenNode;
	while (cur != NULL)
	{
		if (AREEQUAL(cur->name,"center"))
		{
			int segments, n, m, i;
			DOUBLE x, y, xaxis, yaxis, rotate, c, s;
			data_segment *edges;
			data_point *points;
			
			number++;
			if (number > 1)
	    	{
				fprintf(stderr,"Only one center is allowed in regions of type ellipse\n");
				exit(EXIT_FAILURE); 
			}
			n = r->npoints;
			m = r->nsegments;
			x = get_double_value(cur,"x",&found,XML_VALUE_REQUIRED);     
			y = get_double_value(cur,"y",&found,XML_VALUE_REQUIRED);
			xaxis = get_double_value(cur,"ax",&found,XML_VALUE_REQUIRED); 
			yaxis = get_double_value(cur,"ay",&found,XML_VALUE_REQUIRED); 
			segments = get_int_value(cur,"segments",&found,XML_VALUE_REQUIRED);
			if (segments < 3)
				segments = 3;
			rotate = get_double_value(cur,"rotate",&found,XML_VALUE_OPTIONAL);
			if (found == 0)
				rotate = 0.0;
			make_vector(edges,segments);
			make_vector(points,segments);
			c = cos(PI*rotate/180);
			s = sin(PI*rotate/180);
			for(i = 0;i < segments;i++)
			{
				DOUBLE t = 2*i*PI/segments;
				DOUBLE p, q;
				points[i].point_no = n + i;
				p = xaxis * cos(t);
				q = yaxis * sin(t);	
				points[i].x = x + c * p - s * q;
				points[i].y = y + s * p + c * q;
				points[i].bc = bc;
				edges[i].segment_no = m + i;
				edges[i].point_no_1 = n + i;
				edges[i].point_no_2 = n + (i + 1) % segments;
				edges[i].bc = bc;
			}
			append_to_vector(r->points,points,n,segments,FREE_APPENDED);
			append_to_vector(r->segments,edges,m,segments,FREE_APPENDED);
			r->npoints = segments + n;
			r->nsegments = segments + m;
		}
		cur = cur->next;
	}
	if(number == 0)
	{
		fprintf(stderr,"One center is required in regions of type ellipse\n");
		exit(EXIT_FAILURE); 
	}
	return;

error_malloc:
	fprintf(stderr,"Error in a malloc or calloc call\n");
	exit(EXIT_FAILURE);	
}

static DOUBLE **parseBezier(xmlNodePtr node,int *npoints)
{
	DOUBLE **data;
	xmlNodePtr cur;
	int alloc, i, found;

	alloc = 256;
	i = 0;
	make_matrix(data,2,alloc);
	cur =  node->xmlChildrenNode;
	while (cur != NULL)
    {
		if (AREEQUAL(cur->name,"point"))
		{
			DOUBLE x, y;
			x = get_double_value(cur,"x",&found,XML_VALUE_REQUIRED);    
			y = get_double_value(cur,"y",&found,XML_VALUE_REQUIRED);
			if (i == alloc)
			{
				alloc += 256;
				expand_vector(data[0],alloc);
				expand_vector(data[1],alloc);
			}
			data[0][i] = x;
			data[1][i] = y;
			i++;
		}
		cur = cur->next;
	}
	*npoints = i;
	return data;

error_malloc:
	fprintf(stderr,"Error in a malloc or calloc call\n");
	exit(EXIT_FAILURE);	
}

DOUBLE *bezier_curve(DOUBLE **data,int begin,int end,DOUBLE t)
{
	DOUBLE *p;
	make_vector(p,2);

	if (begin == end)
	{
		p[0] = data[0][begin];
		p[1] = data[1][begin];
		return p;
	}
	DOUBLE *r, *s;
	r = bezier_curve(data,begin,end - 1,t);
	s = bezier_curve(data,begin + 1,end,t);

	p[0] = (1 - t) * r[0] + t * s[0];
	p[1] = (1 - t) * r[1] + t * s[1];; 
	free_vector(r);
	free_vector(s);
	return p;

error_malloc:
	fprintf(stderr,"Error in a malloc or calloc call\n");
	exit(EXIT_FAILURE);	
}


static void parsePath(xmlNodePtr node,DataRegion r)
{
	xmlChar *hole;
	xmlNodePtr cur;
	int found, n, m, size, alloc, i, pbc;
	data_segment *edges;
	data_point *points;

	hole =  get_string_value(node,"ishole",&found,XML_VALUE_OPTIONAL);
	if (found && AREEQUAL(hole,"true"))
	{
		data_hole *h;
		DOUBLE x, y;
		int n;
		n = r->nholes;
		x = get_double_value(node,"x",&found,XML_VALUE_REQUIRED);
		y = get_double_value(node,"y",&found,XML_VALUE_REQUIRED);
		make_vector(h,1);
		h->x = x;
		h->y = y;
		append_to_vector(r->holes,h,n,1,FREE_APPENDED);
		r->nholes = n + 1;
	}
	
	pbc = get_int_value(node,"bc",&found,XML_VALUE_OPTIONAL);
	if (found == 0)
		pbc = FEM_BC_DIRICHLET;

	n = r->npoints;
	m = r->nsegments;
	alloc = 256;
	make_vector(points,alloc);
	i = 0;
	
	cur =  node->xmlChildrenNode;
	while (cur != NULL)
	{
		if (AREEQUAL(cur->name,"point"))
		{
			DOUBLE x, y;
			int bc, sbc;
			x = get_double_value(cur,"x",&found,XML_VALUE_REQUIRED);    
			y = get_double_value(cur,"y",&found,XML_VALUE_REQUIRED);
			bc = get_int_value(cur,"bc",&found,XML_VALUE_OPTIONAL);
			if (found == 0)
				bc = pbc;
			sbc = get_int_value(cur,"sbc",&found,XML_VALUE_OPTIONAL);
			if (found == 0)
				sbc = bc;

			if (i == alloc)
	    	{
	      		alloc += 256;
	      		expand_vector(points,alloc);
	    	}
			points[i].point_no = n + i;
			points[i].x = x;
			points[i].y = y;
			points[i].bc = bc;
			points[i].sbc = sbc;
			i++;
		}
		else if(AREEQUAL(cur->name,"bezier"))
		{
			DOUBLE **bezier, *p;
			int nbezier, segments, k, bc, sbc;
			
			segments = get_int_value(cur,"segments",&found,XML_VALUE_REQUIRED);
			if (segments < 5)
				segments = 5;
			bc = get_int_value(cur,"bc",&found,XML_VALUE_OPTIONAL);
			if (found == 0)
				bc = pbc;
			sbc = get_int_value(cur,"sbc",&found,XML_VALUE_OPTIONAL);
			if (found == 0)
				sbc = bc;			
			bezier = parseBezier(cur,&nbezier);
			for (k = 1;k <= segments;k++)
	    	{
	      		DOUBLE t;
	      		t = (DOUBLE)(1.0 * k) / segments;
	      		p = bezier_curve(bezier,0,nbezier - 1,t);
	      		if (i == alloc)
				{
					alloc += 256;
					expand_vector(points,alloc);
				}
				points[i].point_no = n + i;
				points[i].x = p[0];
				points[i].y = p[1];
				points[i].bc = bc;
				points[i].sbc = sbc;

				i++;
				free_vector(p);
			}
			free_matrix(bezier);
		}
		else if(AREEQUAL(cur->name,"ellipticarc"))
		{
			int segments, k;
			DOUBLE x0, y0, ax, ay;
			DOUBLE rotate, begin, end, delta, t, bc, sbc;
			DOUBLE cr, sr;
			segments = get_int_value(cur,"segments",&found,XML_VALUE_REQUIRED);
			if (segments < 5)
				segments = 5;
			bc = get_int_value(cur,"bc",&found,XML_VALUE_OPTIONAL);
			if (found == 0)
				bc = pbc;
			sbc = get_int_value(cur,"sbc",&found,XML_VALUE_OPTIONAL);
			if (found == 0)
				sbc = bc;		
			x0 = get_double_value(cur,"x",&found,XML_VALUE_REQUIRED);
			y0 = get_double_value(cur,"y",&found,XML_VALUE_REQUIRED);
			ax = get_double_value(cur,"ax",&found,XML_VALUE_REQUIRED);
			ay = get_double_value(cur,"ay",&found,XML_VALUE_REQUIRED);
			rotate = get_double_value(cur,"rotate",&found,XML_VALUE_OPTIONAL);
			if (found == 0)
				rotate = 0.0;
			begin = get_double_value(cur,"begin",&found,XML_VALUE_REQUIRED);
			end = get_double_value(cur,"end",&found,XML_VALUE_REQUIRED);
			rotate *= PI/180.0;
			begin *= PI/180.0;
			end *= PI/180.0;
			delta = (end - begin) / segments;
			cr = cos(rotate);
			sr = sin(rotate);
			for (k=1;k <= segments;k++)
			{
				DOUBLE x, y, x1, y1;
				t = begin + k * delta;
				if (i == alloc)
				{
					alloc += 256;
					expand_vector(points,alloc);
				}
				x1 = ax * cos(t);
				y1 = ay * sin(t);
				x = cr * x1 - sr * y1;
				y = sr * x1 + cr * y1;
				x += x0;
				y += y0;
				points[i].point_no = n + i;
				points[i].x = x;
				points[i].y = y;
				points[i].bc = bc;
				points[i].sbc = sbc;
				i++;
			}
		}
		else if(AREEQUAL(cur->name,"curve"))
		{
			int segments, k;
			DOUBLE begin, end, delta, bc, sbc, t;
			xmlChar *valuex, *valuey;
			void *fx, *fy;
			char *names[] = { "t" };
			DOUBLE valuet[1];
			segments = get_int_value(cur,"segments",&found,XML_VALUE_REQUIRED);
			if (segments < 5)
				segments = 5;
			bc = get_int_value(cur,"bc",&found,XML_VALUE_OPTIONAL);
			if (found == 0)
				bc = pbc;
			sbc = get_int_value(cur,"sbc",&found,XML_VALUE_OPTIONAL);
			if (found == 0)
				sbc = bc;		
			begin = get_double_value(cur,"begin",&found,XML_VALUE_REQUIRED);
			end = get_double_value(cur,"end",&found,XML_VALUE_REQUIRED);
			valuex = get_string_value(cur,"x",&found,XML_VALUE_REQUIRED);
			valuey = get_string_value(cur,"y",&found,XML_VALUE_REQUIRED);
			fx = evaluator_create(valuex);
			fy = evaluator_create(valuey);
			if ((fx == NULL) || (fy == NULL))
			{
				fprintf(stderr,"Invalid evaluator\n");
				exit(EXIT_FAILURE);
			}
			delta = (end - begin) / segments;
			for (k=1;k <= segments;k++)
			{
				DOUBLE x, y;
				t = begin + k * delta;
				if (i == alloc)
				{
					alloc += 256;
					expand_vector(points,alloc);
				}
				valuet[0] = t;
				x = evaluator_evaluate(fx,1,names,valuet);
				y = evaluator_evaluate(fy,1,names,valuet);
				points[i].point_no = n + i;
				points[i].x = x;
				points[i].y = y;
				points[i].bc = bc;
				points[i].sbc = sbc;
				i++;
			}
			evaluator_destroy(fx);
			evaluator_destroy(fy);
			xmlFree(valuex);
			xmlFree(valuey);
		}
		cur = cur->next;
	}
	/*
		Append edges
	*/
	size = i;
	r->npoints = n + size;  
 	make_vector(edges,size);
	
	for (i = 0;i < size;i++)
	{
		edges[i].segment_no = m + i;
		edges[i].point_no_1 = n + i;
		edges[i].point_no_2 = n + ((i + 1) % size);
		edges[i].bc = points[i].sbc;
	}
	append_to_vector(r->points,points,n,size,FREE_APPENDED);
	append_to_vector(r->segments,edges,m,size,FREE_APPENDED);
	r->nsegments = m + size;
	return;

error_malloc:
	fprintf(stderr,"Error in a malloc or calloc call\n");
	exit(EXIT_FAILURE);
}

static void parseRegion(xmlNodePtr node,DataRegion r)
{
	xmlChar *type;
	int found;
	type = get_string_value(node,"type",&found,XML_VALUE_REQUIRED);
	if(AREEQUAL(type,"polygon"))
		parsePolygon(node,r);
	else if(AREEQUAL(type,"circle"))
		parseCircle(node,r);
	else if(AREEQUAL(type,"square"))
		parseSquare(node,r);
	else if(AREEQUAL(type,"ellipse"))
		parseEllipse(node,r);
	else if(AREEQUAL(type,"path"))
		parsePath(node,r);
  
	xmlFree(type);
}

static void parseHoles(xmlNodePtr node,DataRegion r)
{
	xmlNodePtr cur;
	cur =  node->xmlChildrenNode;
	while (cur != NULL)
	{
		if (AREEQUAL(cur->name,"hole"))
		{
			data_hole *h;
			DOUBLE x, y;
			int n, found;
			n = r->nholes;
			x = get_double_value(cur,"x",&found,XML_VALUE_REQUIRED);
			y = get_double_value(cur,"y",&found,XML_VALUE_REQUIRED);
			make_vector(h,1);
			h->x = x;
			h->y = y;
			append_to_vector(r->holes,h,n,1,FREE_APPENDED);
			r->nholes = n+1;
		}
		cur = cur->next;
	}
	return;
	
error_malloc:
	fprintf(stderr,"Error in a malloc or calloc call\n");
	exit(EXIT_FAILURE);
}

DataRegion parseXMLRegionDocument(char *docname)
{
	xmlDocPtr doc;
	xmlNodePtr cur;
 	DataRegion r;
 	int error;
	
	doc = NULL;
	cur = NULL;
	r = NULL;
	error = 1;
  
	if ((r = data_region_new()) == NULL)
		return NULL;
	r->npoints = 0;
	r->nsegments = 0;
	r->nholes = 0;
	r->points = NULL;
	
	r->segments = NULL;
	r->holes = NULL;
  
	if ((doc = xmlParseFile(docname)) == NULL)
		goto final;
	
	if ((cur = xmlDocGetRootElement(doc)) == NULL)
		goto final;
	
	if (! AREEQUAL(cur->name,"regions"))
		goto final;
	
	cur = cur->xmlChildrenNode;
	while (cur != NULL)
	{
		if(AREEQUAL(cur->name,"region"))
			parseRegion(cur,r);
		if(AREEQUAL(cur->name,"holes"))
			parseHoles(cur,r);
		cur = cur->next;
	}
	error = 0;
	
final:
	if (error == 1)
	{
		freeDataRegion(r);
		r = NULL;
	}
  	if (doc != NULL)
 		xmlFreeDoc(doc);
	return r;
}

void region_to_asy(DataRegion r,const char *outfile)
{
	FILE *fp;
	time_t now;
	int i;
	
	if ((fp = fopen(outfile, "w")) == NULL)
	{
		fprintf(stderr, "cannot open file %s for writing\n", outfile);
		return;
	}

	fprintf(fp, "/*\n");
	fprintf(fp, " * Title: (%s)\n", outfile);
	fprintf(fp, " * Creator: Study and implementation in C of the Finite Element Method in dimension 2, %s\n", __FILE__);
	now = time(NULL);
	fprintf(fp, " * CreationDate: %s", ctime(&now));
	fprintf(fp,"*/\n");
	fprintf(fp,"unitsize(1.0cm);\n");

	for(i = 0;i < r->nsegments;i++)
	{
		int begin, end;
		begin = r->segments[i].point_no_1;
		end = r->segments[i].point_no_2;
		if (r->segments[i].bc == FEM_BC_NEUMANN)
		{
			fprintf(fp,"draw((%g,%g)--(%g,%g),red+linewidth(1pt));\n",
					r->points[begin].x,
					r->points[begin].y,
					r->points[end].x,
					r->points[end].y);
		}
		else if (r->segments[i].bc == FEM_BC_ROBIN)
		{
			fprintf(fp,"draw((%g,%g)--(%g,%g),green+linewidth(1pt));\n",
					r->points[begin].x,
					r->points[begin].y,
					r->points[end].x,
					r->points[end].y);
		}
		else
		{
			fprintf(fp,"draw((%g,%g)--(%g,%g),blue+linewidth(1pt));\n",
					r->points[begin].x,
					r->points[begin].y,
					r->points[end].x,
					r->points[end].y);
		}
	}
	
	for (i = 0;i < r->nholes;i++)
	{
		fprintf(fp,"dot((%g,%g),gray+linewidth(5pt));\n",
				r->holes[i].x,
				r->holes[i].y); 
	}
	fclose(fp);
}

void print_region(DataRegion r)
{
	int i;
	printf("Number of points: %d\n",r->npoints);
	printf("Number of segmens: %d\n",r->nsegments);
	printf("Number of holes: %d\n",r->nholes);
	
	printf("\nPoints:\n");
	for (i = 0;i < r->npoints;i++)
		printf("  Point %d   (%.12g,%.12g)\n",i,r->points[i].x,r->points[i].y);
	
	printf("\nSegments:\n");
	for (i = 0;i < r->nsegments;i++)
		printf("  %d--%d\n",r->segments[i].point_no_1,r->segments[i].point_no_2);
	
	printf("\nHoles:\n");
	for (i = 0;i < r->nholes;i++)
		printf("  (%.12g,%.12g)\n",r->holes[i].x,r->holes[i].y);
}

