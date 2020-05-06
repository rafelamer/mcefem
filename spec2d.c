/**************************************************************************************
* Filename:   spec2d.c
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

#define AREEQUAL(a,b)   !xmlStrcmp((a),(const xmlChar *)(b))
#define LCAREEQUAL(a,b)   !xmlStrcasecmp((a),(const xmlChar *)(b))

void free_spec2D(Specification2D *spec)
{
	if (*spec == NULL)
		goto final;
	freeLagrange1DValues((*spec)->lv1d);
	freeLagrange2DValues((*spec)->lv2d,(*spec)->degree);
	freeDataMesh((*spec)->mesh);
	freeFunctionTable((*spec)->funs);
	if ((*spec)->Lua != NULL)
		lua_close((*spec)->Lua);
	if ((*spec)->pModule != NULL)
		Py_Finalize_Functions((*spec)->pModule);

final:
	free(*spec);
	*spec = NULL;
}

static int parseRegion(xmlNodePtr node,Specification2D spec,DataRegion *rg,DOUBLE *size)
{
	int found, error; 
	xmlChar *value;

  	error = 1;
  	value = NULL;
	value = get_string_value(node,"filename",&found,XML_VALUE_REQUIRED);
	if ((*rg = parseXMLRegionDocument((char *)value)) == NULL)
		goto final;

	*size = get_double_value(node,"size",&found,XML_VALUE_REQUIRED);
	
	error = 0;

final:
	if (value != NULL)
		xmlFree(value);
	return error;
}

static int parseSettings(xmlNodePtr node,Specification2D spec)
{
	xmlNodePtr cur;
	int l = 0, q1d = 0, q2d = 0, i = 0;

	cur = node->xmlChildrenNode;
	while (cur != NULL)
	{
		int found;
		if (AREEQUAL(cur->name,"lagrange"))
		{
			spec->degree = get_int_value(cur,"degree",&found,XML_VALUE_REQUIRED);
			if (spec->degree < 1)
				spec->degree = 1;
			l++;
		}
		if (AREEQUAL(cur->name,"quadrature1d"))
		{
			spec->npoints1d = get_int_value(cur,"npoints",&found,XML_VALUE_REQUIRED);
			if (spec->npoints1d < 3)
				spec->npoints1d = 3;
			q1d++;
		}
		if (AREEQUAL(cur->name,"quadrature2d"))
		{
			spec->npoints2d = get_int_value(cur,"npoints",&found,XML_VALUE_REQUIRED);
			if (spec->npoints2d < 3)
				spec->npoints2d = 3;
			q2d++;
		}


		if ((l > 1) || (q1d > 1) || (q2d > 1))
			return 1;
		cur = cur->next;
	}
	if ((l != 1) || (q1d != 1) || (q2d != 1))
		return 1;
	spec->qdat1d = gauss_data(&(spec->npoints1d));
	if ((spec->lv1d = LagrangeAtGaussPoints(spec->degree,spec->qdat1d,spec->npoints1d)) == NULL)
		return 1;
	int degree = 1;
	spec->qdat2d = get_twb_data(&degree,&(spec->npoints2d));
	if ((spec->lv2d = LagrangeAtTWBPoints(spec->degree,spec->qdat2d,spec->npoints2d)) == NULL)
		return 1;

	return 0;
}

static char *strstrip(char *s)
{
	size_t size;
	char *end, *begin;

	begin = s;
	size = strlen(s);

	if (size == 0)
	return  begin;

	end = begin + size - 1;
	while (end >= begin && isspace(*end))
		end--;
	*(end + 1) = '\0';

	while (*begin && isspace(*begin))
		begin++;

	return begin;
}

static int parseFunctions(xmlDocPtr doc,xmlNodePtr node,Specification2D spec)
{
	xmlNodePtr cur;
	int a = 0, b1 = 0, b2 = 0, c = 0, f = 0;
	int d = 0, n = 0, rA = 0, rB = 0, rC = 0; 

	cur = node->xmlChildrenNode;
	while (cur != NULL)
	{
		int error = 0;
		if (AREEQUAL(cur->name,"text"))
		{
			cur = cur->next;
			continue;
		}
		if(AREEQUAL(cur->name,"python"))
		{
			xmlChar *name;
			int found;
			name = get_string_value(cur,"filename",&found,XML_VALUE_OPTIONAL);
			spec->pModule = Py_Initialize_Functions("tfgfem",name);
			xmlFree(name);
			if (spec->pModule == NULL)
				return 0;
		}
		if (AREEQUAL(cur->name,"function"))
		{
			xmlChar *type, *name;
			int found;
			type = get_string_value(cur,"type",&found,XML_VALUE_REQUIRED);
			name = get_string_value(cur,"name",&found,XML_VALUE_REQUIRED);
			if (LCAREEQUAL(type,"lua"))
			{
				xmlChar *value;
				char *aux;
				if (spec->Lua == NULL)
				{
					spec->Lua = luaL_newstate();
					luaL_openlibs(spec->Lua);
				}
				value = get_string_content(doc,cur);
				aux = strstrip((char *)value);
				add_luafunction(&(spec->funs),name,spec->Lua,aux);
				xmlFree(value);
			}
			else if (LCAREEQUAL(type,"c"))
			{
				FunctionTable item;
				HASH_FIND_STR(spec->funs,name,item);
				if (item == NULL)
					error = 1;
				if (item->type != CFUNCTION2D)
					error = 1;
			}
			else if (LCAREEQUAL(type,"matheval"))
			{
				xmlChar *value;
				char *aux;
				value = get_string_content(doc,cur);
				aux = strstrip((char *)value);
				add_mathevaluator(&(spec->funs),name,aux);
				xmlFree(value);
			}
			else if (LCAREEQUAL(type,"python"))
			{
				add_pythonfunction(&(spec->funs),name,spec->pModule);
			}
			if (AREEQUAL(name,"a"))
				a++;
			else if (AREEQUAL(name,"b1"))
				b1++;
			else if (AREEQUAL(name,"b2"))
				b2++;
			else if (AREEQUAL(name,"c"))
				c++;
			else if (AREEQUAL(name,"f"))
				f++;
			else if (AREEQUAL(name,"d"))
				d++;
			else if (AREEQUAL(name,"n"))
				n++;
			else if (AREEQUAL(name,"A"))
				rA++;
			else if (AREEQUAL(name,"B"))
				rB++;
			else if (AREEQUAL(name,"C"))
				rC++;			
			else
				error = 1;

			xmlFree(type);
			xmlFree(name);
			if (error == 1)
				return 1;
		}
		if ((a > 1) || (b1 > 1) || (b2 > 1) || (c > 1) || (f > 1) ||
			(d > 1) || (n > 1) || (rA > 1) || (rB > 1) || (rC > 1))
			return 1;
		cur = cur->next;
	}

	if ((a != 1) || (b1 != 1) || (b2 != 1) || (c != 1) || (f != 1))
		return 1;

	return 0;
}

static void parseOutput(xmlNodePtr node,Specification2D spec)
{
	int found;
	xmlChar *value; 

	value = get_string_value(node,"type",&found,XML_VALUE_REQUIRED);
	strncpy(spec->type,value,63);
	xmlFree(value);
	value = get_string_value(node,"filename",&found,XML_VALUE_OPTIONAL);
	strncpy(spec->filename,value,1023);
	xmlFree(value);
	spec->elementvalues = get_int_value(node,"elementvalues",&found,XML_VALUE_REQUIRED);
}


Specification2D parseXMLSpec2DDocument(char *docname,FunctionTable funs)
{
	xmlDocPtr doc;
	xmlNodePtr cur;
	Specification2D spec;
	int error = 1;
	int r = 0, s = 0, f = 0, o = 0;

	if ((spec = (Specification2D)malloc(sizeof(spec2D))) == NULL)
		goto final;
	spec->funs = funs;
	spec->pModule = NULL;
	spec->Lua = NULL;
	DataRegion rg;
	DOUBLE size;

	if ((doc = xmlParseFile(docname)) == NULL)
		goto final;

	if ((cur = xmlDocGetRootElement(doc)) == NULL)
		goto final;

	if (! AREEQUAL(cur->name,"fem2d"))
		goto final;

	cur = cur->xmlChildrenNode;
	while (cur != NULL)
	{
		if (AREEQUAL(cur->name,"region"))
		{
			r++;
			if ((error = parseRegion(cur,spec,&rg,&size)) == 1)
				goto final;
		}
		if (AREEQUAL(cur->name,"settings"))
		{
			s++;
			if ((error = parseSettings(cur,spec)) == 1)
				goto final;
		}
		if (AREEQUAL(cur->name,"functions"))
		{
			f++;
			if ((error = parseFunctions(doc,cur,spec)) == 1)
				goto final;
		}
		if (AREEQUAL(cur->name,"output"))
		{
			o++;
			parseOutput(cur,spec);
		}

		if ((r > 1) || (s > 1) || (f > 1) || (o > 1))
			goto final;
		cur = cur->next;
	}
	if ((r != 1) || (s != 1) || (f != 1))
		goto final;

	DataMesh mesh;
	if ((mesh = make_mesh(rg,size)) == NULL)
		goto final;
	spec->mesh = mesh;

	FunctionTable item;
	if (spec->mesh->hasDirichletBC)
	{
		HASH_FIND_STR(spec->funs,"d",item);
		if (item == NULL)
			goto final;
	}
	if (spec->mesh->hasNeumannBC)
	{
		HASH_FIND_STR(spec->funs,"n",item);
		if (item == NULL)
			goto final;
	}
	if (spec->mesh->hasRobinBC)
	{
		HASH_FIND_STR(spec->funs,"A",item);
		if (item == NULL)
			goto final;
		HASH_FIND_STR(spec->funs,"B",item);
		if (item == NULL)
			goto final;
		HASH_FIND_STR(spec->funs,"C",item);
		if (item == NULL)
			goto final;
	}

	error = 0;

final:
	if (rg != NULL)
		freeDataRegion(rg);
	if ((error == 1) && (spec != NULL))
		freeSpecification2D(spec);
	if (doc != NULL)
		xmlFreeDoc(doc);
	return spec;
}
