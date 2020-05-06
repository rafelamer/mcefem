/**************************************************************************************
* Filename:   xml.c
* Author:     
* Copyright:  
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
#include <unistd.h>

#define GETVALUE(n,a) xmlGetProp((n),(const xmlChar *)a)

/*
  Functions for reading attrinutes from XML files
 */
int get_int_value(xmlNodePtr node,char *name,int *found,int optional)
{
	xmlChar *str;
	char *ptr;
	int value;
	str = GETVALUE(node,name);
	if ((optional != 0) && (str == NULL))
	{
		fprintf(stderr,"Value %s required but not set in %s\n",name,node->name);
		exit(EXIT_FAILURE);
	}
	if (str == NULL)
	{
		*found = 0;
		return 0;
	}
	*found = 1;
	value = strtol((const char *)str,&ptr,10);
	if (*ptr != '\0')
	{
		fprintf(stderr,"Incorrect value of %s: %s\n",name,str);
		exit(EXIT_FAILURE);
	}
	xmlFree(str);
	return value;
}


DOUBLE get_double_value(xmlNodePtr node,char *name,int *found,int optional)
{
	xmlChar *str;
	char *ptr;
	DOUBLE value;
	str = GETVALUE(node,name);
	if ((optional != 0) && (str == NULL))
	{
		fprintf(stderr,"Value %s required but not set in %s\n",name,node->name);
		exit(EXIT_FAILURE);
	}
	if (str == NULL)
	{
		*found = 0;
		return 0.0;
	}
	*found = 1;
	value = strtod((const char *)str,&ptr);
	if (*ptr != '\0')
	{
		fprintf(stderr,"Incorrect value of %s: %s\n",name,str);
		exit(EXIT_FAILURE);
	}
	xmlFree(str);
	return value;
}

xmlChar *get_string_value(xmlNodePtr node,char *name,int *found,int optional)
{
	xmlChar *str;
	str = GETVALUE(node,name);
	if ((optional != 0) && (str == NULL))
	{
		fprintf(stderr,"Value %s required but not set in %s\n",name,node->name);
		exit(EXIT_FAILURE);
	}
	*found = 1;
	if (str == NULL)
		*found = 0;
	return str;
}

xmlChar *get_string_content(xmlDocPtr doc,xmlNodePtr node)
{
	xmlChar *str;
	str = xmlNodeListGetString(doc,node->xmlChildrenNode,1);
	return str;
}