import axes;
import graph;

usepackage("times");
usepackage("mtpro2");

unitsize(1cm);

picture e = axes(0,15,300,320);
picture g = grid(0,15,300,320);
add(g);
add(e);

file in=input("result06.txt").line();
real[][] data=in.dimension(0,0);
data=transpose(data);
draw(graph(data[0],data[1]),blue+linewidth(0.3mm));