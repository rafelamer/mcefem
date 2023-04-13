<<<<<<< HEAD
CFLAGS = -fPIC -O3 -I. -I/usr/include/libxml2 -I/usr/include/python3.10 \
-I/usr/include/petsc
=======
CFLAGS = -fPIC -O3 -I. -I/usr/include/libxml2 -I/usr/include/lua5.2 -I/usr/include/python2.7 \
-I/usr/lib/petscdir/3.12/include/ -I/usr/lib/x86_64-linux-gnu/openmpi/include/
>>>>>>> 5a269dde681c97b5e07b5a4f8af074db92140fcb
CC = gcc
TARGET = libtfgfem.so.1.0.0
OBJECTS = triangle.o sparsem.o gauss.o xml.o regions.o mesh.o \
	twb.o quadrature.o triangulation.o lagrange.o fem1d.o fem2d.o \
	functions.o spec1d.o spec2d.o gsl.o petsc.o

INCLUDES = tfgfem.h

$(TARGET): $(OBJECTS)
	$(CC) -shared -fPIC -Wl,-soname,libtfgfem.so.1 -o $(TARGET) $(OBJECTS)

install: $(TARGET)
	cp $(INCLUDES) /usr/local/include/
	cp $(TARGET) /usr/local/lib/
	ln -sf       /usr/local/lib/$(TARGET) /usr/local/lib/libtfgfem.so.1
	ln -sf       /usr/local/lib/$(TARGET) /usr/local/lib/libtfgfem.so

clean:
	rm -f $(OBJECTS) $(TARGET)
