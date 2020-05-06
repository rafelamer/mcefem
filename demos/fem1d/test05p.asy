import axes;
import graph;

usepackage("times");
usepackage("mtpro2");

unitsize(1.75cm);

picture e = axes(0,10,-3,3);
picture g = grid(0,10,-3,3);
add(g);
add(e);

file in=input("test05.txt").line();
real[][] data=in.dimension(0,0);
data=transpose(data);
draw(graph(data[0],data[1]),red+linewidth(0.3mm));

file in=input("result05p.txt").line();
real[][] data=in.dimension(0,0);
data=transpose(data);
draw(graph(data[0],data[1]),blue+linewidth(0.1mm));