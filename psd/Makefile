CC=gcc
CXX=g++
CFLAGS=-Wall -O3 -Itipsylib
LDFLAGS=-Ltipsylib
LIBS=-ltipsy -lm

psd: psd.cpp
	$(CXX) $(CFLAGS) $(LDFLAGS) psd.cpp $(LIBS) -o psd

clean:
	rm -f psd
	$(MAKE) -C tipsylib clean
