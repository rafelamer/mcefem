#!/bin/bash

for i in 1 2 3 4 5 5p; do
	tfgfem --infile=test0$i.xml --dimension=1
	asy -f pdf test0$i.asy
done
