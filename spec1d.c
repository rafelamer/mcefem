/**************************************************************************************
* Filename:   spec1d.c
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

void free_spec1D(Specification1D *spec)
{
	if (*spec == NULL)
		goto final;
	freeLagrange1DValues((*spec)->lv);
	freeFunctionTable((*spec)->funs);
	if ((*spec)->Lua != NULL)
		lua_close((*spec)->Lua);
	if ((*spec)->pModule != NULL)
		Py_Finalize_Functions((*spec)->pModule);

final:
	free(*spec);
	*spec = NULL;
}

static int parseBoundaries(xmlNodePtr node,Specification1D spec)
{
	xmlNodePtr cur;
	int left = 0, right = 0;
	int error = 1;

	cur = node->xmlChildrenNode;
	while (cur != NULL)
	{
		xmlChar *type;
		int found;
		int index;
		type = NULL;

		if (AREEQUAL(cur->name,"text"))
		{
			cur = cur->next;
			continue;
		}

		if (AREEQUAL(cur->name,"left"))
		{
			type =  get_string_value(cur,"type",&found,XML_VALUE_REQUIRED);
			index = 0;
			left++;
		}
		else if (AREEQUAL(cur->name,"right"))
		{
			type = get_string_value(cur,"type",&found,XML_VALUE_REQUIRED);
			index = 1;
			right++;
		}
		else
			goto final;

		if (LCAREEQUAL(type,"dirichlet"))
		{
			spec->bc[index][0] = get_double_value(cur,"value",&found,XML_VALUE_REQUIRED);
			spec->bctype[index] = FEM_BC_DIRICHLET;
		}
		if (LCAREEQUAL(type,"neumann"))
		{
			spec->bc[index][0] = get_double_value(cur,"value",&found,XML_VALUE_REQUIRED);
			spec->bctype[index] = FEM_BC_NEUMANN;
		}
		if (LCAREEQUAL(type,"robin"))
		{
			spec->bc[index][0] = get_double_value(cur,"A",&found,XML_VALUE_REQUIRED);
			spec->bc[index][1] = get_double_value(cur,"B",&found,XML_VALUE_REQUIRED);
			spec->bc[index][2] = get_double_value(cur,"C",&found,XML_VALUE_REQUIRED);
			spec->bctype[index] = FEM_BC_ROBIN;
		}
		xmlFree(type);
		if ((left > 1) || (right > 1))
			goto final;
		cur = cur->next;
	}
	if ((left != 1) || (right != 1))
		goto final;
	error = 0;
final:
	return error;
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

static int parseFunctions(xmlDocPtr doc,xmlNodePtr node,Specification1D spec)
{
	xmlNodePtr cur;
	int a2 = 0, a1 = 0, a0 = 0, f = 0;

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
				if (item->type != CFUNCTION1D)
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
			if (AREEQUAL(name,"a0"))
				a0++;
			else if (AREEQUAL(name,"a1"))
				a1++;
			else if (AREEQUAL(name,"a2"))
				a2++;
			else if (AREEQUAL(name,"f"))
				f++;
			else
				error = 1;

			xmlFree(type);
			xmlFree(name);
			if (error == 1)
				return 1;
		}
		if ((a2 > 1) || (a1 > 1) || (a0 > 1) || (f > 1))
			return 1;
		cur = cur->next;
	}

	if ((a2 != 1) || (a1 != 1) || (a0 != 1) || (f != 1))
		return 1;
	return 0;
}

static int parseSettings(xmlNodePtr node,Specification1D spec)
{
	xmlNodePtr cur;
	int e = 0, l = 0, q = 0, i = 0;

	cur = node->xmlChildrenNode;
	while (cur != NULL)
	{
		int found;
		if (AREEQUAL(cur->name,"elements"))
		{
			spec->elements = get_int_value(cur,"number",&found,XML_VALUE_REQUIRED);
			if (spec->elements < 1)
				spec->elements = 1;
			e++;
		}
		if (AREEQUAL(cur->name,"lagrange"))
		{
			spec->degree = get_int_value(cur,"degree",&found,XML_VALUE_REQUIRED);
			if (spec->degree < 1)
				spec->degree = 1;
			l++;
		}
		if (AREEQUAL(cur->name,"quadrature"))
		{
			spec->npoints = get_int_value(cur,"npoints",&found,XML_VALUE_REQUIRED);
			if (spec->npoints < 3)
				spec->npoints = 3;
			q++;
		}
		if (AREEQUAL(cur->name,"interval"))
		{
			spec->xa = get_double_value(cur,"xa",&found,XML_VALUE_REQUIRED);
			spec->xb = get_double_value(cur,"xb",&found,XML_VALUE_REQUIRED);
			if (spec->xa >= spec->xb)
				return 1;
			i++;
		}
		if ((e > 1) || (l > 1) || (q > 1) || (i > 1))
			return 1;
		cur = cur->next;
	}
	if ((e != 1) || (l != 1) || (q != 1) || (i != 1))
		return 1;
	spec->qdat = gauss_data(&(spec->npoints));
	if ((spec->lv = LagrangeAtGaussPoints(spec->degree,spec->qdat,spec->npoints)) == NULL)
		return 1;

	return 0;
}

static void parseOutput(xmlNodePtr node,Specification1D spec)
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


Specification1D parseXMLSpec1DDocument(char *docname,FunctionTable funs)
{
	xmlDocPtr doc;
	xmlNodePtr cur;
	Specification1D spec;
	int error = 1;
	int b = 0, f = 0, s = 0, o = 0;

	if ((spec = (Specification1D)malloc(sizeof(spec1D))) == NULL)
		goto final;
	spec->funs = funs;
	spec->pModule = NULL;
	spec->Lua = NULL;

	if ((doc = xmlParseFile(docname)) == NULL)
		goto final;

	if ((cur = xmlDocGetRootElement(doc)) == NULL)
		goto final;

	if (! AREEQUAL(cur->name,"fem1d"))
		goto final;

	cur = cur->xmlChildrenNode;
	while (cur != NULL)
	{
		if (AREEQUAL(cur->name,"boundaries"))
		{
			b++;
			if ((error = parseBoundaries(cur,spec)) == 1)
				goto final;
		}
		if (AREEQUAL(cur->name,"functions"))
		{
			f++;
			if ((error = parseFunctions(doc,cur,spec)) == 1)
				goto final;
		}
		if (AREEQUAL(cur->name,"settings"))
		{
			s++;
			if ((error = parseSettings(cur,spec)) == 1)
				goto final;
		}
		if (AREEQUAL(cur->name,"output"))
		{
			o++;
			parseOutput(cur,spec);
		}

		if ((b > 1) || (f > 1) || (s > 1) || (o > 1))
			goto final;
		cur = cur->next;
	}
	if ((b != 1) || (s != 1) || (f != 1))
		goto final;
	error = 0;

final:
	if ((error == 1) && (spec != NULL))
		freeSpecification1D(spec);
	if (doc != NULL)
		xmlFreeDoc(doc);
	return spec;
}
