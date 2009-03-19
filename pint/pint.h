#ifndef PINT_H
#define PINT_H

#ifndef PINT_MAIN
#define EXTERN extern
#else
#define EXTERN
#endif

#ifndef USE_TREE_JS
#define USE_TREE_JS  0
#endif
#ifndef USE_TREE_JS
#define USE_TREE_JPC 0
#endif

#include <stdio.h>

#include "particle.h"

struct tree;

struct env
{
    struct particle *ps;
    struct particle **p;
    int N;

    struct tree *tree;
};

/*============================================================================
 * Constants
 *==========================================================================*/

#define MAX_PIC 8

#define SAHA_TIME 1

enum 
{ 
    OUT_ENERGY = 0,
    OUT_DENSITY,
    OUT_MASS,
    OUT_ACCEL_HIST,
    OUT_TIMESTEPS,
    OUT_PARTICLE_INFO, // must be last in list
    OUT_NUM
};

enum
{
    SAVE_NOTHING = 0,
    SAVE_TIPSY,
    SAVE_ASCII,
    SAVE_ASCII_ONE_FILE,
};

/*============================================================================
 * Physical Constants
 *==========================================================================*/
#define G (6.6742e-11)                  /* gravitational constant */
//#define K2 (2.959122082855911025e-4)
#define K2 1

/* default softening */
#define DEFAULT_SOFTENING (0)
/* initial time step */
#define TS (3e-6)

/*============================================================================
 * Macros
 *==========================================================================*/

#define VL(__vl) if (verbosity >= __vl)
#define LOG(__ll) if (loglevel >= __ll && logfp != NULL)

#define X_0(exp) do { const int _=0; exp; } while(0)
#define X_1(exp) do { const int _=1; exp; } while(0)
#define X_2(exp) do { const int _=2; exp; } while(0)
#define X3(exp)  do { X_0(exp); X_1(exp); X_2(exp); } while(0)

#define SWAP(A,B,T) { T = A; A = B; B = T; }
//#define SWAP(a,b,T) (((a) == (b)) || (((a) ^= (b)), ((b) ^= (a)), ((a) ^= (b))))
#define PARTITION(P,T,ELEM,i,j,CMPL,CMPU) \
do {\
    while (i <= j && ((P[i] ELEM) CMPL)) { ++i; } \
    while (i <= j && ((P[j] ELEM) CMPU)) { --j; } \
    while (i < j) { \
        SWAP(P[i], P[j], T); \
        while ((P[++i] ELEM) CMPL) { } \
        while ((P[--j] ELEM) CMPU) { } \
    }\
} while(0)

#define SWAP2(A,B,T) { T = A; A = B; B = T; }
#define PARTITION2(A,ELEM,p,r,i,CMPL) \
do {\
    while (i <= j && ((P[i] ELEM) CMPL)) { ++i; } \
    while (i <= j && ((P[j] ELEM) CMPU)) { --j; } \
    while (i < j) { \
        SWAP(P[i], P[j], T); \
        while ((P[++i] ELEM) CMPL) { } \
        while ((P[--j] ELEM) CMPU) { } \
    }\
} while(0)


/* absolute value because abs() has problems */
#define ABS(x) ((x) < 0 ? (-(x)) : (x))
/* the linear interperlated position of a particle at time */
//#define POSAT(p, i, t) ((p)->r[i] - (((p)->time-(t))*(p)->v[i]))
#define POSAT(p, i, t) ((p)->r[i])
#define INTERP_POSITION POSAT

#if 1
#define ACCEL_HIST(a) \
{ \
    int __bin = (int)MMIN(ACCEL_HIST_SIZE-1, log(MAG(a))*10+(ACCEL_HIST_SIZE/2)); \
    accel_hist[(int)MMAX(0, (double)__bin)]++; \
}
#endif
    //fprintf(stderr, "%22.15e %22.15e %i\n", MAG(a), log(MAG(a)), __bin); 

#ifdef ACCEL_HIST
#define ACCEL_HIST_SIZE 400
EXTERN int accel_hist[ACCEL_HIST_SIZE];
#define RESET_ACCEL_HIST() { memset(accel_hist, 0, ACCEL_HIST_SIZE*sizeof(int)); }
#else
#define RESET_ACCEL_HIST()
#endif

/*============================================================================
 * Globals
 *==========================================================================*/
EXTERN int NPARTICLES;

EXTERN struct particle **ps;

/*============================================================================
 * I/O Globals
 *==========================================================================*/
EXTERN FILE *out, *err, *in, *logfp;
EXTERN char *infile;
EXTERN char *outfilebase;
EXTERN unsigned int verbosity;
EXTERN unsigned int loglevel;

/*============================================================================
 * Prototypes
 *==========================================================================*/
real MMAX(real a, real b);
real MMIN(real a, real b);
real SQDIST(struct particle *p, struct particle *q);
real DIST(struct particle *p, struct particle *q);
real SQMAG(real *a);
real MAG(real *a);

#endif 
