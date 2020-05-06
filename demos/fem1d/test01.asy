import axes;
import graph;

usepackage("times");
usepackage("mtpro2");

unitsize(1.75cm);

picture e = axes(0,6,-5,3);
picture g = grid(0,6,-5,3);
add(g);
add(e);

file in=input("test01.txt").line();
real[][] data=in.dimension(0,0);
data=transpose(data);
draw(graph(data[0],data[1]),red+linewidth(0.4mm));

file in=input("result01.txt").line();
real[][] data=in.dimension(0,0);
data=transpose(data);
draw(graph(data[0],data[1]),blue+linewidth(0.1mm));