CC?=gcc
CXX?=g++
WARN=-Wfloat-equal -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Wmissing-noreturn -Winline
INCLUDE=-I$(HOME)/zTools/include/tipsylib
#LIBS=-L$(HOME)/zTools/lib64
#LIBS=-L$(HOME)/zTools/lib
#CPPFLAGS=-g -O3 -Wall $(WARN) $(LISTING) $(INCLUDE) $(LIBS)
CFLAGS=-O3 -Wall $(WARN) $(LISTING) $(INCLUDE) $(LIBS) -g -ftree-vectorize -ftree-vectorizer-verbose=5
CXXFLAGS=$(CFLAGS)


tipsplit: tipsplit.cpp
	$(CXX) $(CXXFLAGS) -o $@ tipsplit.cpp $(HOME)/zTools/lib64/libtipsy.a

tipsubsampl: tipsubsampl.cpp
	$(CXX) $(CXXFLAGS) -o $@ tipsubsampl.cpp $(HOME)/zTools/lib64/libtipsy.a

denminmax: denminmax.cpp
	$(CXX) $(CXXFLAGS) -lpthread -o $@ denminmax.cpp $(HOME)/zTools/lib64/libtipsy.a

denminmax2: denminmax2.cpp
	$(CXX) $(CXXFLAGS) -o $@ denminmax2.cpp $(HOME)/zTools/lib64/libtipsy.a

cutbox: cutbox.cpp
	$(CXX) $(CXXFLAGS) -o $@ cutbox.cpp $(HOME)/zTools/lib64/libtipsy.a

markbox: markbox.cpp
	$(CXX) $(CXXFLAGS) -o $@ markbox.cpp $(HOME)/zTools/lib64/libtipsy.a

mkptgrid: mkptgrid.cpp
	$(CXX) $(CXXFLAGS) -o $@ mkptgrid.cpp $(HOME)/zTools/lib64/libtipsy.a

numden: numden.cpp
	$(CXX) $(CXXFLAGS) -o $@ numden.cpp $(HOME)/zTools/lib64/libtipsy.a

minmaxpot: minmaxpot.cpp
	$(CXX) $(CXXFLAGS) -o $@ minmaxpot.cpp $(HOME)/zTools/lib64/libtipsy.a

tipsytimes: tipsytimes.cpp
	$(CXX) $(CXXFLAGS) -o $@ tipsytimes.cpp $(HOME)/zTools/lib64/libtipsy.a

swapposvel: swapposvel.cpp
	$(CXX) $(CXXFLAGS) -o $@ swapposvel.cpp $(HOME)/zTools/lib64/libtipsy.a

tipexplode: tipexplode.cpp
	$(CXX) $(CXXFLAGS) -o $@ tipexplode.cpp $(HOME)/zTools/lib64/libtipsy.a

ahf2tipgrp: ahf2tipgrp.c
	$(CC) $(CFLAGS) -o $@ ahf2tipgrp.c

ahf2grp: ahf2grp.c
	$(CC) $(CFLAGS) -o $@ ahf2grp.c

all: tipsplit tipsubsampl denminmax denminmax2 cutbox mkptgrid numden minmaxpot ahf2tipgrp 

