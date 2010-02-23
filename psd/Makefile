CC=gcc
CXX=g++
CFLAGS=-Wall -O2 -Itipsylib
LDFLAGS=-Ltipsylib
LIBS=-ltipsy -lm
TIPSYLIB=libtipsy.a

psd: psd.cpp $(TIPSYLIB)
	$(CXX) $(CFLAGS) $(LDFLAGS) psd.cpp $(LIBS) -o psd

$(TIPSYLIB):
	$(MAKE) -C tipsylib

clean:
	rm -f psd
	$(MAKE) -C tipsylib clean
