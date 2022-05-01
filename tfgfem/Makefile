CFLAGS = -O3 -I. -I/usr/include/libxml2 -I/usr/lib/petscdir/3.12/include/ -I/usr/lib/x86_64-linux-gnu/openmpi/include/ $(shell python2.7-config --cflags) $(shell python2.7-config --ldflags)
LIBS = -ltfgfem -lumfpack -lm -lmatheval -lxml2 -llua5.2 -lut -lpython2.7 -lgsl -lpetsc -lmpi
LDFLAGS = -L/usr/local/lib/
CC = gcc
TARGET = tfgfem
OBJECTS = tfgfem.o cmdline.o
INCLUDES = cmdline.h

$(TARGET): $(OBJECTS) $(INCLUDES) tfgfem.c
	$(CC) -o $(TARGET) $(OBJECTS) $(LDFLAGS) $(LIBS)

$(OBJECTS): tfgfem.c cmdline.c

cmdline.c: cmdline.ggo
	gengetopt --input=cmdline.ggo

all: $(TARGET)

install: $(TARGET)
	cp $(TARGET) /usr/local/bin/

clean:
	rm -f $(OBJECTS) cmdline.[ch] $(TARGET) 
