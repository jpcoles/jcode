%module pytipsy
%{
#include "ftipsy.hpp"
#include "tipsypos.h"
%}

%include <stdint.i>

%typemap(memberin) float [ANY] {
  memcpy($1, $input, $1_dim0 * sizeof(float));
}

%typemap(out) float [ANY] {
  int i;
  $result = PyList_New($1_dim0);
  for (i = 0; i < $1_dim0; i++) {
    PyObject *o = PyFloat_FromDouble((double) $1[i]);
    PyList_SetItem($result,i,o);
  }
}

%typemap(out) tipsypos::section_type {
  $result = PyInt_FromLong($1);
}

%rename(read)  iTipsy::operator>>;
%rename(write) oTipsy::operator<<;

%feature("shadow") iTipsy::operator>> %{
    def __iter__(self):
        return self
    def next(self):
        p = self.read()
        if p is None: raise StopIteration()
        return p

    def read(self, *args):
        if len(args) > 0: return $action(self, *args)
        pos = self.tellg().section()
        if pos == tipsypos.gas:
            g = TipsyGasParticle()
            self.read(g)
            return g
        if pos == tipsypos.dark:
            d = TipsyDarkParticle()
            $action(self, d)
            #self.read(d)
            return d
        if pos == tipsypos.star:
            s = TipsyStarParticle()
            self.read(s)
            return s
        if pos == tipsypos.header:
            h = TipsyHeader()
            self.read(h)
            return h
        return None
%}

/*
%extend iTipsy {
    TipsyBaseParticle *read() 
    { 
        switch($self->tellg().section())
        {
            case tipsypos::header:
            {
                TipsyHeader *h = new TipsyHeader();
                *$self >> *h;
                return (TipsyBaseParticle *)h;
                break;
            }
            case tipsypos::gas:
            {
                TipsyGasParticle *g = new TipsyGasParticle();
                *$self >> *g;
                return (TipsyBaseParticle *)g;
                break;
            }
            case tipsypos::dark:
            {
                TipsyDarkParticle *d = new TipsyDarkParticle();
                *$self >> *d;
                return (TipsyBaseParticle *)d;
                break;
            }
            case tipsypos::star:
            {
                TipsyStarParticle *s = new TipsyStarParticle();
                *$self >> *s;
                return (TipsyBaseParticle *)s;
                break;
            }
        }

        return NULL;
    }
};
*/

%include "ftipsy.hpp"

%rename(__equal__) tipsypos::operator==;
%rename(__assign__) tipsypos::operator=;
%rename(section_ptr) tipsypos::section();
%rename(offset_ptr) tipsypos::offset();
%extend tipsypos {
    section_type section(void) { return $self->section(); }
    offset_type  offset(void)  { return $self->offset();  }
}

%include "tipsypos.h"


