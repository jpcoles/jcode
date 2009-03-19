CXX=g++
CXXFLAGS=-g -Wall -I$(HOME)/zTools/include/tipsylib
LDFLAGS=-Ltipsylib -ltipsy -ljpeg

tipmri: tipmri.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) tipmri.cpp tipsylib/libtipsy.a -o tipmri

clean:
	rm -f tipmri
	$(MAKE) -C tipsylib clean
