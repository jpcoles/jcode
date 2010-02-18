#include <iostream>
#include <stdlib.h>
#include <stdint.h>
#include <float.h>      // FLT_MIN, FLT_MAX
#include <math.h>
#include "timing.h"
#include "ftipsy.hpp"

CPUDEFS

using namespace std;

//typedef unsigned int uint128_t;
typedef unsigned int uint128_t __attribute__((mode(TI)));

typedef struct
{
    float x[3], v[3], M;
    float rho;
} particle_t;

typedef struct
{
    uint128_t key;
    size_t pid;
} Key;

int _key_compar(const void *a, const void *b)
{
    return ((Key *)a)->key < ((Key *)b)->key;
}

int main(int argc, char **argv)
{
    ifTipsy           in;
    TipsyHeader        h;
    TipsyGasParticle   g;
    TipsyDarkParticle  d;
    TipsyStarParticle  s;
    uint32_t i, j;
    double t0,t1;

    float dx=0, dv=0;

    cerr << "sizeof(uint128_t)=" << sizeof(uint128_t) << endl;

    if (argc < 2)
    {
        cerr << "Usage: psd <tipsy-file> [dx [dv]]" << endl;
        exit(2);
    }

    if (argc >= 3) dx = atof(argv[2]);
    if (argc >= 4) dv = atof(argv[3]);

    in.open(argv[1], "standard");
    if (!in.is_open())
    {
        cerr << "Can't open " << argv[1] << endl;
        exit(1);
    }

    in >> h;
    uint32_t N = h.h_nBodies;
    float    a = h.h_time;

    Key *keys = (Key *)malloc(N * sizeof(Key));
    if (keys == NULL)
    {
        in.close();
        cerr << "Not enough memory for keys." << endl;
        exit(1);
    }

    particle_t *ps = (particle_t *)malloc(N * sizeof(particle_t));
    if (ps == NULL)
    {
        in.close();
        cerr << "Not enough memory for particles." << endl;
        exit(1);
    }

#define M_SQRT8PI_3 sqrt(8 * M_PI / 3);

    //dx *= a;
    //dv *= a;

    //--------------------------------------------------------------------------
    cerr << "Loading snapshot ("
         << "N=" << N
         << "; mem=" << (N*(sizeof(particle_t) + sizeof(Key))) << " bytes"
         << "; a=" << h.h_time
         << ")... ";
    //--------------------------------------------------------------------------
#define READ(_p, _n)                    \
do {                                    \
    double _t, _t0 = CPUTIME;           \
    uint32_t last_i=0;                  \
    double timeout = 5;                 \
    _t = _t0;                           \
    for (i=0; i < _n; i++, j++) {       \
        in >> _p;                       \
        ps[j].x[0] = _p.pos[0] * a;     \
        ps[j].x[1] = _p.pos[1] * a;     \
        ps[j].x[2] = _p.pos[2] * a;     \
        ps[j].v[0] = _p.vel[0] * a;     \
        ps[j].v[1] = _p.vel[1] * a;     \
        ps[j].v[2] = _p.vel[2] * a;     \
        ps[j].M    = _p.mass;     \
        if ((CPUTIME)-_t > timeout) {    \
            _t = CPUTIME;                \
            if (!timed_out) cerr << endl; \
            timed_out = 1;              \
            cerr << "\r" << ((j*100.0)/N) << "%" \
                 << " (" << (i/(_t-_t0)) << "/s)                  ";\
            last_i = i; \
            timeout = 2;                \
            }} \
    if ((CPUTIME)-_t > timeout) {         \
        timed_out = 1;              \
        if (!timed_out) cerr << endl; \
            cerr << "\r" << ((j*100.0)/N) << "%" \
                 << " (" << (i/(_t-_t0)) << "/s)                  ";\
        } \
} while (0)

    j=0;
    t0 = CPUTIME;
    int timed_out=0;
    READ(g, h.h_nSph);
    READ(d, h.h_nDark);
    READ(s, h.h_nStar);
    in.close();
    t1 = CPUTIME;
    if (timed_out) cerr << endl << "Load time: ";

    cerr << (t1-t0) << "s" << endl;

    //--------------------------------------------------------------------------
    cerr << endl << "Pass 1 (Find range)... "; cerr.flush();
    //--------------------------------------------------------------------------
    double start_time = CPUTIME;

    float xrange[6] = {FLT_MAX,FLT_MAX,FLT_MAX, FLT_MIN,FLT_MIN,FLT_MIN}, 
          vrange[6] = {FLT_MAX,FLT_MAX,FLT_MAX, FLT_MIN,FLT_MIN,FLT_MIN};

    t0 = start_time;
    for (i=0; i < N; i++)
    {
        if (ps[i].x[0] < xrange[0]) xrange[0] = ps[i].x[0];
        if (ps[i].x[0] > xrange[3]) xrange[3] = ps[i].x[0];
        if (ps[i].x[1] < xrange[1]) xrange[1] = ps[i].x[1];
        if (ps[i].x[1] > xrange[4]) xrange[4] = ps[i].x[1];
        if (ps[i].x[2] < xrange[2]) xrange[2] = ps[i].x[2];
        if (ps[i].x[2] > xrange[5]) xrange[5] = ps[i].x[2];

        if (ps[i].v[0] < vrange[0]) vrange[0] = ps[i].v[0];
        if (ps[i].v[0] > vrange[3]) vrange[3] = ps[i].v[0];
        if (ps[i].v[1] < vrange[1]) vrange[1] = ps[i].v[1];
        if (ps[i].v[1] > vrange[4]) vrange[4] = ps[i].v[1];
        if (ps[i].v[2] < vrange[2]) vrange[2] = ps[i].v[2];
        if (ps[i].v[2] > vrange[5]) vrange[5] = ps[i].v[2];
    }
    t1 = CPUTIME;
    cerr << (t1-t0) << "s" << endl;


    //--------------------------------------------------------------------------
    // Find the extents of the simulation in position and velocity space.
    //--------------------------------------------------------------------------

    float max_xrange = max( xrange[3]-xrange[0], 
                       max( xrange[4]-xrange[1],
                            xrange[5]-xrange[2] ));

    float max_vrange = max( vrange[3]-vrange[0], 
                       max( vrange[4]-vrange[1],
                            vrange[5]-vrange[2] ));

    if (max_xrange == 0) 
    {
        cerr << "Empty range over all position dimensions" << endl;
        exit(1);
    }

    if (max_vrange == 0) 
    {
        cerr << "Empty range over all velocity dimensions" << endl;
        exit(1);
    }

    //--------------------------------------------------------------------------
    // Compute the number of cells in each dimension if values of dx and dv are
    // given. If they were not given, calculate them assuming 1000 cells in 
    // each dimension.
    //--------------------------------------------------------------------------
    //uint32_t nxcells = 1000;
    //uint32_t nvcells = 1000;

    uint32_t nxcells = pow(h.h_nBodies, 1./3);
    uint32_t nvcells = pow(h.h_nBodies, 1./3);

    if (dx != 0) nxcells = (uint32_t)ceil(max_xrange * 1.001 / dx);
    if (dv != 0) nvcells = (uint32_t)ceil(max_vrange * 1.001 / dv);

    if (dx == 0) dx = max_xrange * 1.0001 / nxcells;
    if (dv == 0) dv = max_vrange * 1.0001 / nvcells;

    //--------------------------------------------------------------------------
    cerr << "Pass 2 (Generate keys;"
         << " dx=" << dx 
         << " dv=" << dv 
         << " nxcells=" << nxcells
         << " nvcells=" << nvcells
         << ")... "; cerr.flush();
    //--------------------------------------------------------------------------

    uint128_t ncells1 = nxcells;
    uint128_t ncells2 = ncells1 * nxcells;
    uint128_t ncells3 = ncells2 * nxcells;
    uint128_t ncells4 = ncells3 * nvcells;
    uint128_t ncells5 = ncells4 * nvcells;

    t0 = CPUTIME;
    for (i=0; i < N; i++)
    {
        keys[i].key = ((uint128_t)((ps[i].x[0] - xrange[0]) / dx))
                    + ((uint128_t)((ps[i].x[1] - xrange[1]) / dx)) * ncells1
                    + ((uint128_t)((ps[i].x[2] - xrange[2]) / dx)) * ncells2
                    + ((uint128_t)((ps[i].v[0] - vrange[0]) / dv)) * ncells3
                    + ((uint128_t)((ps[i].v[1] - vrange[1]) / dv)) * ncells4
                    + ((uint128_t)((ps[i].v[2] - vrange[2]) / dv)) * ncells5
                    ;
        keys[i].pid = i;
    }
    t1 = CPUTIME;
    cerr << (t1-t0) << "s" << endl;
    
    //--------------------------------------------------------------------------
    cerr << "Pass 3 (Sort)... "; cerr.flush();
    //--------------------------------------------------------------------------
    t0 = CPUTIME;
    qsort(keys, N, sizeof(Key), _key_compar);
    t1 = CPUTIME;
    cerr << (t1-t0) << "s" << endl;


    //--------------------------------------------------------------------------
    cerr << "Pass 4 (Integrate)... "; cerr.flush();
    //--------------------------------------------------------------------------
    t0 = CPUTIME;

    double    S=0, S0=0;
    double    rho;
    uint32_t  longest_run=0, shortest_run=INT_MAX;
    uint32_t  f=0;
    uint32_t  nkeys=0;
    uint32_t  first_in_run=0;
    uint128_t key = -1;
    for (i=0; i < N; i++)
    {
        if (key != keys[i].key)
        {
            if (f != 0)
            {
                double f0 = ((double)f) / (double)N;
                S  += f0 * log(f0);
                S0 += f0;

                if (f > longest_run)  longest_run  = f;
                if (f < shortest_run) shortest_run = f;

                rho = 0;
                for (j=first_in_run; j < i; j++)
                    rho += ps[keys[j].pid].M;

                rho /= pow(dx,3) * pow(dv,3);

                for (j=first_in_run; j < i; j++)
                {
                    ps[keys[j].pid].rho = rho;
                }
            }



            f = 0;
            key = keys[i].key;
            nkeys++;

            first_in_run = i;

        }

        f++;
    }

    if (f != 0)
    {
        double f0 = ((double)f) / (double)N;
        S  += f0 * log(f0);
        S0 += f0;
        if (f > longest_run)  longest_run  = f;
        if (f < shortest_run) shortest_run = f;
    }

    t1 = CPUTIME;
    cerr << (t1-t0) << "s" << endl;

    //--------------------------------------------------------------------------

    double end_time = CPUTIME;
    float r = (double)nkeys / N;

    cerr << "Stats:" << endl;
    cerr << "\tCompute time:      " << (end_time - start_time) << "s" << endl
         << "\tncells:            " << nkeys << endl
         << "\tncells/N:          " << r << endl
         << "\t<ps/cell>:         " << (1/r) << endl
         << "\tRuns (long,short): " << longest_run << "," << shortest_run << endl
         << "\tUnweighted sum:    " << S0 << endl
         ;

    //S = -S / (pow(dx,3)*pow(dv,3));
    S = -S;
    //cout << "S = " << S << endl;
    cerr << "S = " << S << endl;

    cout << h.h_nBodies << endl;
    for (i=0; i < N; i++)
    {
        cout << ps[i].rho << endl;
    }

    return 0;
}

