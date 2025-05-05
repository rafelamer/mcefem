CFLAGS = -O3 -I. -I/usr/include/libxml2 -I/usr/include/petsc $(shell python3-config --cflags) $(shell python3-config --ldflags)
LIBS = -ltfgfem -lumfpack -lm -lmatheval -lxml2 -llua -lpython3.13 -lgsl -lpetsc
LDFLAGS = -L/usr/local/lib 
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
