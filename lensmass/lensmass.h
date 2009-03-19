#ifndef LENSMASS_H
#define LENSMASS_H

#ifdef __MACOSX__
#include <Accelerate/Accelerate.h>
typedef __CLPK_integer integer;
typedef __CLPK_real real;
#else
#ifdef __LINUX__
#undef complex
#define Skip_f2c_Undefs
#include "f2c.h"
extern "C"
{
#include "clapack.h"
}
#define __CLPK_integer integer
#define __CLPK_real real
#endif
#endif

#define SHOW_MATRIX_FLAT        0x01
#define SHOW_MATRIX             0x02
#define SHOW_EIGEN_VALUES       0x04
#define SHOW_EIGEN_VECTORS_R    0x20
#define SHOW_EIGEN_VECTORS_L    0x40

#define SHOW_DEFAULT                \
(SHOW_MATRIX_FLAT | SHOW_EIGEN_VALUES | SHOW_EIGEN_VECTORS_R)

typedef struct eigen_s {
    double re;
    double im;
    float *lvec;
    float *rvec;
    int dim;
} eigen_t;

#endif

