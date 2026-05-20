import axes;
import graph;

usepackage("times");
usepackage("mtpro2");

unitsize(1cm);

picture e = axes(-10,10,-10,12);
picture g = grid(-10,10,-10,12);
add(g);
add(e);

file in=input("test07.txt").line();
real[][] data=in.dimension(0,0);
data=transpose(data);
draw(graph(data[0],data[1]),red+linewidth(0.3mm));

file in=input("result07.txt").line();
real[][] data=in.dimension(0,0);
data=transpose(data);
draw(graph(data[0],data[1]),blue+linewidth(0.15mm));