CC=gcc
UNAME=$(shell uname)
USERSRC=$(HOME)/src

#CFLAGS=-g -O -pg -lm -DNDEBUG 
CFLAGS=-g -std=c99 -O3 -lm 
#CFLAGS=-Wall -lm -g -O3 -pg 
PINT_SRC_C=pint.c io.c tree_js.c tree_jpc.c
PINT_SRC_H=particle.h io.h pint.h tree_js.h tree_jpc.h
OPTIONS=-I$(USERSRC)/glfw-2.5/include

ifeq ($(CC), pgcc)
	CFLAGS += -B
endif
ifeq ($(CC), cc)
	#CFLAGS += -O0 -pg -Wall -std=c99 -fprofile-arcs -ftest-coverage
	CFLAGS += -Wall -std=c99
endif

ifeq ($(UNAME), Darwin)
        OPTIONS += -L$(USERSRC)/glfw-2.5/lib/macosx
        LDOPTIONS += -lglfw -framework OpenGL -framework AGL -framework Carbon

        #OPTIONS += -I$(USERSRC)/OpenGL/Headers/
endif

ifeq ($(UNAME), Linux)
        #LIBS += -lglut
endif

jpc: $(PINT_SRC_C) $(PINT_SRC_H) tipsy.o Makefile
	$(CC) $(CFLAGS) $(OPTIONS) $(LDOPTIONS) $(LIBS) -DUSE_TREE_JPC=1 $(PINT_SRC_C) tipsy.o -o pint 

js: $(PINT_SRC_C) $(PINT_SRC_H) tipsy.o Makefile
	$(CC) $(CFLAGS) $(OPTIONS) $(LDOPTIONS) $(LIBS) -DUSE_TREE_JS=1 $(PINT_SRC_C) tipsy.o -o pint 

tipsy.o: tipsy.c tipsy.h 
	$(CC) $(CFLAGS) -c tipsy.c

timestep: timestep.c Makefile
	$(CC) $(CFLAGS) timestep.c -o timestep

clean:
	rm -f pint *.o

.PHONY: jpc js
