#!/usr/bin/python3
"""This demo program solves Poisson's equation

    - div grad u(x, y) = f(x, y)

on the unit square with source f given by

    f(x, y) = 10*exp(-((x - 0.5)^2 + (y - 0.5)^2) / 0.02)

and boundary conditions given by

    u(x, y) = 0        for x = 0 or x = 1
du/dn(x, y) = sin(5*x) for y = 0 or y = 1
"""

from dolfin import *

mesh = UnitSquareMesh(48, 48)
V = FunctionSpace(mesh, "Lagrange", 2)

def boundary(x):
    return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS

u0 = Constant(0.0)
bc = DirichletBC(V, u0, boundary)

u = TrialFunction(V)
v = TestFunction(V)
f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)",degree=2)
g = Expression("sin(5*x[0])",degree=2)
a = inner(grad(u), grad(v))*dx
L = f*v*dx + g*v*ds

u = Function(V)
solve(a == L, u, bc)

values = u.compute_vertex_values()
count = 0
for i in mesh.coordinates():
	X = i[0]
	Y = i[1]
	Z = values[count]
	count += 1
	print (X, Y, Z)

