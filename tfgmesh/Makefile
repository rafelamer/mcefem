CFLAGS = -O3 -I. -I/usr/include/libxml2 -I/usr/lib/petscdir/3.7.3/x86_64-linux-gnu-real/include/ -I/usr/lib/openmpi/include/ $(shell python-config --cflags) $(shell python-config --ldflags)
LIBS = -ltfgfem -lumfpack -lm -lmatheval -lxml2 -llua5.2 -lut -lpython2.7 -lgsl -lpetsc -lmpi
LDFLAGS = -L/usr/local/lib/
CC = gcc
TARGET = tfgmesh
OBJECTS = tfgmesh.o cmdline.o
INCLUDES = cmdline.h

$(TARGET): $(OBJECTS) $(INCLUDES) tfgmesh.c
	$(CC) -o $(TARGET) $(OBJECTS) $(LDFLAGS) $(LIBS)

$(OBJECTS): tfgmesh.c cmdline.c

cmdline.c: cmdline.ggo
	gengetopt --input=cmdline.ggo

all: $(TARGET)

install: $(TARGET)
	cp $(TARGET) /usr/local/bin/

clean:
	rm -f $(OBJECTS) cmdline.[ch] $(TARGET)
