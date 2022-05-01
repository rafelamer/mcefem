/**************************************************************************************
* Filename:   functions.c
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

PyObject *Py_Initialize_Functions(char *programname,char *filename)
{
	PyObject *pModule, *pName;

	Py_SetProgramName(programname);
	Py_Initialize();
	PySys_SetPath(".");
	pName = PyString_FromString(filename);
	pModule = PyImport_Import(pName);

	Py_DECREF(pName);
	return pModule;
}

void Py_Finalize_Functions(PyObject *module)
{
	Py_DECREF(module);
	Py_Finalize();
}

void add_cfunction1d(FunctionTable *hash,char *name,FunctionX f)
{
	FunctionTable new;
	int error = 1;

	if ((new = (FunctionTable)malloc(sizeof(function_table))) == NULL)
		goto final;

	strncpy(new->name,name,127);
	new->type = CFUNCTION1D;
	new->f1d = f;
	HASH_ADD_STR(*hash,name,new);
	error = 0;

final:
	if (error == 1)
	{
		if (new != NULL)
			free(new);
		fprintf(stderr,"Error generating a C function %s\n",name);
		exit(EXIT_FAILURE);
	}
}

void add_cfunction2d(FunctionTable *hash,char *name,FunctionXY f)
{
	FunctionTable new;
	int error = 1;

	if ((new = (FunctionTable)malloc(sizeof(function_table))) == NULL)
		goto final;

	strncpy(new->name,name,127);
	new->type = CFUNCTION2D;
	new->f2d = f;
	HASH_ADD_STR(*hash,name,new);
	error = 0;

final:
	if (error == 1)
	{
		if (new != NULL)
			free(new);
		fprintf(stderr,"Error generating a C function %s\n",name);
		exit(EXIT_FAILURE);
	}
}

void add_mathevaluator(FunctionTable *hash,char *name,char *function)
{
	FunctionTable new;
	int error = 1;

	if ((new = (FunctionTable)malloc(sizeof(function_table))) == NULL)
		goto final;

	strncpy(new->name,name,127);
	new->type = MATHEVALUATOR;
	new->evaluator = evaluator_create(function);
	if (new->evaluator == NULL)
		goto final;
	HASH_ADD_STR(*hash,name,new);
	error = 0;

final:
	if (error == 1)
	{
		if (new != NULL)
			free(new);
		fprintf(stderr,"Error generating a math evaluator for function %s\n",function);
		exit(EXIT_FAILURE);
	}
}

void add_luafunction(FunctionTable *hash,char *name,lua_State *Lua,char *function)
{
	FunctionTable new;
	int error = 1;

	if ((new = (FunctionTable)malloc(sizeof(function_table))) == NULL)
		goto final;

	strncpy(new->name,name,127);
	new->type = LUAFUNCTION;
	new->Lua = Lua;

	if (luaL_loadbuffer(Lua,function,strlen(function),name) != 0)
		goto final;
	if (lua_pcall(Lua, 0, 0, 0) != 0)
		goto final;
	HASH_ADD_STR(*hash,name,new);
	error = 0;

final:
	if (error == 1)
	{
		if (new != NULL)
			free(new);
		fprintf(stderr,"Error generating the Lua function %s\n%s\n",name,function);
		exit(EXIT_FAILURE);
	}
}

void add_pythonfunction(FunctionTable *hash,char *name,PyObject *python)
{
	FunctionTable new;
	PyObject *pFunc;
	int error = 1;

	if ((new = (FunctionTable)malloc(sizeof(function_table))) == NULL)
		goto final;

	strncpy(new->name,name,127);
	new->type = PYTHONFUNCTION;
	pFunc = PyObject_GetAttrString(python,name);
	if ((pFunc == NULL) || (!PyCallable_Check(pFunc)))
		goto final;
	new->Python = pFunc;
	HASH_ADD_STR(*hash,name,new);
	error = 0;

final:
	if (error == 1)
	{
		Py_XDECREF(pFunc);
		Py_DECREF(python);
		Py_Finalize();
		if (new != NULL)
			free(new);
		fprintf(stderr,"Error generating the Python function %s\n",name);
		exit(EXIT_FAILURE);
	}
}

DOUBLE callFunctionTable1D(FunctionTable hash,char *name,DOUBLE x)
{
	FunctionTable item;

	HASH_FIND_STR(hash,name,item);
	if (item == NULL)
		goto final;
	if (item->type == MATHEVALUATOR)
		return evaluator_evaluate_x(item->evaluator,x);
	if (item->type == LUAFUNCTION)
	{
		DOUBLE r;
		lua_getglobal(item->Lua,item->name);
		if (lua_isnil(item->Lua, -1))
			goto final;
		lua_pushnumber(item->Lua,x);
		lua_pcall(item->Lua,1,1,0);
		r = lua_tonumber(item->Lua,-1);
		lua_pop(item->Lua,1);
		return r;
	}
	if (item->type == CFUNCTION1D)
		return (item->f1d)(x);
	{
		PyObject *pArgs, *pValue;
		pArgs = PyTuple_New(1);
		if ((pValue = PyFloat_FromDouble(x)) == NULL)
			goto final;
		PyTuple_SetItem(pArgs,0,pValue);
		pValue = PyObject_CallObject(item->Python,pArgs);
		Py_DECREF(pArgs);

		DOUBLE y = PyFloat_AsDouble(pValue);
		Py_DECREF(pValue);
		return y;
	}

	return NAN;

final:
	fprintf(stderr,"Error calling the table function %s\n",name);
	exit(EXIT_FAILURE);
}

DOUBLE callFunctionTable2D(FunctionTable hash,char *name,DOUBLE x,DOUBLE y)
{
	FunctionTable item;

	HASH_FIND_STR(hash,name,item);
	if (item == NULL)
		goto final;
	if (item->type == MATHEVALUATOR)
		return evaluator_evaluate_x_y(item->evaluator,x,y);
	if (item->type == LUAFUNCTION)
	{
		DOUBLE r;
		lua_getglobal(item->Lua,item->name);
		if (lua_isnil(item->Lua, -1))
			goto final;
		lua_pushnumber(item->Lua,x);
		lua_pushnumber(item->Lua,y);
		lua_pcall(item->Lua,2,1,0);
		r = lua_tonumber(item->Lua,-1);
		lua_pop(item->Lua,1);
		return r;
	}
	if (item->type == CFUNCTION1D)
		return (item->f2d)(x,y);
	{
		PyObject *pArgs, *pValue;
		pArgs = PyTuple_New(2);
		if ((pValue = PyFloat_FromDouble(x)) == NULL)
			goto final;
		PyTuple_SetItem(pArgs,0,pValue);
		if ((pValue = PyFloat_FromDouble(y)) == NULL)
			goto final;
		PyTuple_SetItem(pArgs,1,pValue);
		pValue = PyObject_CallObject(item->Python,pArgs);
		Py_DECREF(pArgs);

		DOUBLE r = PyFloat_AsDouble(pValue);
		Py_DECREF(pValue);
		return r;
	}

	return NAN;

final:
	fprintf(stderr,"Error calling the table function %s\n",name);
	exit(EXIT_FAILURE);
}

DOUBLE runFunctionTable1D(FunctionTable item,DOUBLE x)
{
	if (item->type == MATHEVALUATOR)
		return evaluator_evaluate_x(item->evaluator,x);
	if (item->type == LUAFUNCTION)
	{
		DOUBLE r;
		lua_getglobal(item->Lua,item->name);
		if (lua_isnil(item->Lua, -1))
			goto final;
		lua_pushnumber(item->Lua,x);
		lua_pcall(item->Lua,1,1,0);
		r = lua_tonumber(item->Lua,-1);
		lua_pop(item->Lua,1);
		return r;
	}
	if (item->type == CFUNCTION1D)
		return (item->f1d)(x);
	if (item->type == PYTHONFUNCTION)
	{
		PyObject *pArgs, *pValue;
		pArgs = PyTuple_New(1);
		if ((pValue = PyFloat_FromDouble(x)) == NULL)
			goto final;
		PyTuple_SetItem(pArgs,0,pValue);
		pValue = PyObject_CallObject(item->Python,pArgs);
		Py_DECREF(pArgs);

		DOUBLE y = PyFloat_AsDouble(pValue);
		Py_DECREF(pValue);
		return y;
	}

	return NAN;

final:
	fprintf(stderr,"Error calling the table function %s\n",item->name);
	exit(EXIT_FAILURE);
}

DOUBLE runFunctionTable2D(FunctionTable item,DOUBLE x,DOUBLE y)
{
	if (item->type == MATHEVALUATOR)
		return evaluator_evaluate_x_y(item->evaluator,x,y);
	if (item->type == LUAFUNCTION)
	{
		DOUBLE r;
		lua_getglobal(item->Lua,item->name);
		if (lua_isnil(item->Lua, -1))
			goto final;
		lua_pushnumber(item->Lua,x);
		lua_pushnumber(item->Lua,y);
		lua_pcall(item->Lua,2,1,0);
		r = lua_tonumber(item->Lua,-1);
		lua_pop(item->Lua,1);
		return r;
	}
	if (item->type == CFUNCTION1D)
		return (item->f2d)(x,y);
	if (item->type == PYTHONFUNCTION)
	{
		PyObject *pArgs, *pValue;
		pArgs = PyTuple_New(2);
		if ((pValue = PyFloat_FromDouble(x)) == NULL)
			goto final;
		PyTuple_SetItem(pArgs,0,pValue);
		if ((pValue = PyFloat_FromDouble(y)) == NULL)
			goto final;
		PyTuple_SetItem(pArgs,1,pValue);
		pValue = PyObject_CallObject(item->Python,pArgs);
		Py_DECREF(pArgs);

		DOUBLE r = PyFloat_AsDouble(pValue);
		Py_DECREF(pValue);
		return r;
	}

	return NAN;

final:
	fprintf(stderr,"Error calling the table function %s\n",item->name);
	exit(EXIT_FAILURE);
}

void free_function_table(FunctionTable *table)
{
	FunctionTable item, tmp;

	if(*table == NULL)
		return;

	HASH_ITER(hh,*table,item,tmp)
	{
		if (item->type == MATHEVALUATOR)
			evaluator_destroy(item->evaluator);
		if (item->type == PYTHONFUNCTION)
			Py_XDECREF(item->Python);
		HASH_DEL(*table,item);
		free(item);
	}
	*table = NULL;
}
