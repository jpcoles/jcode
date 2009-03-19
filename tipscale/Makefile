CXX?=g++
WARN=-Wfloat-equal -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Wmissing-noreturn -Winline
INCLUDE=-I$(HOME)/zTools/include/tipsylib
CPPFLAGS=-O3 -Wall $(WARN) $(LISTING) $(INCLUDE) $(LIBS)

tipscale: tipscale.cpp
	$(CXX) $(CPPFLAGS) -o $@ tipscale.cpp $(HOME)/zTools/lib64/libtipsy.a

all: tipscale


