#CC=gcc
CC=pgCC
CFLAGS=-g -Wall -O3
CLIBS=-lreadline
LAPACKLIBS=
PLATFORM=$(shell uname)
EXT=$(value CONDOR_EXT)

ifeq ($(PLATFORM),Darwin)
	CFLAGS+=-faltivec -bind_at_load -D__MACOSX__
	LAPACKLIBS+=-framework Accelerate 
endif
ifeq ($(PLATFORM),Linux)
	CFLAGS+= -D__LINUX__
	CLIBS+=-lm 
	ifeq ($(CC),pgCC)
		CFLAGS+=-B
		LAPACKLIBS+=-fPIC -llapack -lblas -lg2c -lpgftnrtl
	endif
	ifeq ($(CC),gcc)
		LAPACKLIBS+=-fPIC /usr/lib64/liblapack.so.3 -lblas -lg2c
	endif
endif

all: lensmass ev2cntr

lensmass: lensmass.c
	$(CC) $(CFLAGS) $(CLIBS) $(LAPACKLIBS) lensmass.c -o lensmass$(EXT)
	cp lensmass$(EXT) lensmass

lensmass-test: lensmass.c
	$(CC) $(CFLAGS) $(CLIBS) $(LAPACKLIBS) -DTEST lensmass.c -o lensmass$(EXT)
	cp lensmass$(EXT) lensmass

ev2cntr: ev2cntr.c
	$(CC) $(CFLAGS) $(CLIBS) ev2cntr.c -o ev2cntr$(EXT)
	cp ev2cntr$(EXT) ev2cntr

