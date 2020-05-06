import axes;
import graph;

usepackage("times");
usepackage("mtpro2");

unitsize(1.75cm);

picture e = axes(-4,4,0,7);
picture g = grid(-4,4,0,7);
add(g);
add(e);

file in=input("test04.txt").line();
real[][] data=in.dimension(0,0);
data=transpose(data);
draw(graph(data[0],data[1]),red+linewidth(0.4mm));

file in=input("result04.txt").line();
real[][] data=in.dimension(0,0);
data=transpose(data);
draw(graph(data[0],data[1]),blue+linewidth(0.1mm));