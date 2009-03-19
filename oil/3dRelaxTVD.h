#ifndef _3DRELAXTVD_H
#define _3DRELAXTVD_H

#define USE_FLOAT 1
#define USE_DOUBLE 0

#if USE_FLOAT
typedef float_t real;
#define pow  powf
#define max  fmaxf
#define min  fminf
#define sqrt sqrtf
#define abs  fabsf
#endif

#if USE_DOUBLE
typedef double_t real;
#define pow  pow
#define max  fmax
#define min  fmin
#define sqrt sqrt
#define abs  fabs
#endif

typedef struct
{
    real p;     /* rho */
    real pvx;
    real pvy;
    real pvz;
    real e;
} cell_t;

typedef int (*limiter_fn_t)(cell_t *f, cell_t *a, cell_t *b);

#define nc 64
#define hc (nc / 2)

static int ic_sedovtaylor();
static void write_output();
static int timestep();
static int averageflux(cell_t *u, cell_t *w, real *c);
static int relaxing(cell_t *u);
static int sweepx();
static int sweepy();
static int sweepz();
static int vanleer(cell_t *f, cell_t *a, cell_t *b);
static int minmod(cell_t *f, cell_t *a, cell_t *b);
static int superbee(cell_t *f, cell_t *a, cell_t *b);

#endif
