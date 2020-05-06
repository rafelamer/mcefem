import math

def a0(x):
	if x <= 2.0:
		return -0.1*((x-1)**2-0.5)
	if x <= 4.0:
		return -0.1*((x-3)**2-0.5)
	if x <= 6.0:
		return -0.1*((x-5)**2-0.5)
	if x <= 8.0:
		return -0.1*((x-7)**2-0.5)
	if x <= 10.0:
		return -0.1*((x-9)**2-0.5)
	return 0.0

def f(x):
	return 0.5*math.floor(x+1)