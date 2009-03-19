# gnu compiler is recommended as it has fast macros for endian conversion
#
#	make TARGET_ARCH="-m32 -march=i686"
#	make TARGET_ARCH="-m64 -march=opteron"
#

# This greatly benifits g++ which assumes i386 on 32 bit machines
#TARGET_ARCH=-march=i686
#TARGET_ARCH=-march=pentium4
#TARGET_ARCH=-march=prescott

#LISTING=-g -Wa,-adhln=$*.lst

CXX=g++
WARN=-Wall -Wfloat-equal -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Wmissing-noreturn -Winline
CPPFLAGS=-fPIC -O3 -Wall -Wextra $(WARN) $(LISTING)

#CXX=icc
#CPPFLAGS=-Kc++ -O3 -xP -g $(WARN) $(LISTING)

#CXX=icpc
#CPPFLAGS=-O3 $(LISTING)

#CXX=pgCC
#CPPFLAGS=-fast

SWIG=swig
SWIGI=pytipsy.i
SWIGPYTHONINCLUDE=-I/usr/include/python2.4

OBJS=ftipsy.o adapter.o native.o standard.o vtipsy.o binner.o zstream.o grid.o
EXTRAOBJS=tdump.o tdiff.o boxstat.o makemark.o makegrid.o gridunion.o tostd.o vtest.o tbin.o
SRCS=ftipsy.cpp adapter.cpp native.cpp standard.cpp vtipsy.cpp zstream.cpp boxstat.cpp makemark.cpp makegrid.cpp gridunion.cpp grid.cpp tdump.cpp tdiff.cpp tostd.cpp vtest.cpp binner.cpp tbin.cpp
HDRS=ftipsy.hpp adapter.h native.h standard.h tipsypos.h vtipsy.hpp zstream.hpp binner.hpp grid.hpp

libtipsy.a : $(OBJS)
	ar r libtipsy.a $(OBJS)

pytipsy : libtipsy.a
	$(SWIG) -Wall -c++ -python $(SWIGI)
	$(CXX) -fPIC -c pytipsy_wrap.cxx $(SWIGPYTHONINCLUDE)
	$(CXX) -fPIC -shared pytipsy_wrap.o libtipsy.a -o _pytipsy.so

all:	libtipsy.a boxstat makemark makegrid gridunion tdump tdiff tostd vtest tbin

clean		:
	rm -f libtipsy.a $(OBJS) $(EXTRAOBJS) boxstat makemark makegrid gridunion tdump tdiff tostd vtest tbin pytipsy_wrap.cxx pytipsy_wrap.o

dep		:
	makedepend -Y. 2>/dev/null $(SRCS) $(HDRS)

boxstat	: boxstat.o libtipsy.a
	$(CXX) $(CPPFLAGS) -o $@ $< -L. -ltipsy

makemark : makemark.o libtipsy.a
	$(CXX) $(CPPFLAGS) -o $@ $< -L. -ltipsy

makegrid : makegrid.o libtipsy.a
	$(CXX) $(CPPFLAGS) -o $@ $< -L. -ltipsy

gridunion : gridunion.o libtipsy.a
	$(CXX) $(CPPFLAGS) -o $@ $< -L. -ltipsy

tbin	: tbin.o libtipsy.a
	$(CXX) $(CPPFLAGS) -o $@ $< -L. -ltipsy


tdump	: tdump.o libtipsy.a
	$(CXX) $(CPPFLAGS) -o $@ $< -L. -ltipsy

tdiff	: tdiff.o libtipsy.a
	$(CXX) $(CPPFLAGS) -o $@ $< -L. -ltipsy

tostd	: tostd.o libtipsy.a
	$(CXX) $(CPPFLAGS) -o $@ $< -L. -ltipsy

vtest	: vtest.o libtipsy.a
	$(CXX) $(CPPFLAGS) -o $@ $< -L. -ltipsy

# DO NOT DELETE

ftipsy.o: ftipsy.hpp tipsypos.h native.h adapter.h standard.h
adapter.o: adapter.h tipsypos.h
native.o: native.h adapter.h tipsypos.h tipsyrec.h
standard.o: standard.h native.h adapter.h tipsypos.h tipsyrec.h
vtipsy.o: vtipsy.hpp ftipsy.hpp tipsypos.h
zstream.o: zstream.hpp
boxstat.o: ftipsy.hpp tipsypos.h tiplim.hpp grid.hpp
makemark.o: ftipsy.hpp tipsypos.h grid.hpp
makegrid.o: grid.hpp
gridunion.o: grid.hpp
grid.o: grid.hpp tiplim.hpp
tdump.o: regex.hpp ftipsy.hpp tipsypos.h
tdiff.o: ftipsy.hpp tipsypos.h
tostd.o: ftipsy.hpp tipsypos.h
vtest.o: vtipsy.hpp ftipsy.hpp tipsypos.h
binner.o: binner.hpp ftipsy.hpp tipsypos.h
tbin.o: ftipsy.hpp tipsypos.h binner.hpp
ftipsy.o: tipsypos.h
adapter.o: tipsypos.h
native.o: adapter.h tipsypos.h
standard.o: native.h adapter.h tipsypos.h
vtipsy.o: ftipsy.hpp tipsypos.h
binner.o: ftipsy.hpp tipsypos.h
