CXX?=g++
WARN=-Wfloat-equal -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Wmissing-noreturn -Winline
INCLUDE=-I$(HOME)/zTools/include/tipsylib
#LIBS=-L$(HOME)/zTools/lib64
#LIBS=-L$(HOME)/zTools/lib
CPPFLAGS=-O3 -Wall $(WARN) $(LISTING) $(INCLUDE) $(LIBS)
#CPPFLAGS=-O3 -Wall $(WARN) $(LISTING) -g -pg #-ftree-vectorize -ftree-vectorizer-verbose=5


tsort: 
	$(CXX) $(CPPFLAGS) -o $@ tsort.cpp $(HOME)/zTools/lib64/libtipsy.a

all: tsort


