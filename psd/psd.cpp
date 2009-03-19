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
    float x[3], v[3];
} particle_t;

int _key_compar(const void *a, const void *b)
{
    return *((uint128_t *)a) < *((uint128_t *)b);
}

int main(int argc, char **argv)
{
    ifTipsy           in;
    TipsyHeader        h;
    TipsyGasParticle   g;
    TipsyDarkParticle  d;
    TipsyStarParticle  s;
    uint32_t i, j;

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

    uint128_t *keys = (uint128_t *)malloc(N * sizeof(uint128_t));
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

    //--------------------------------------------------------------------------
    cout << "Loading snapshot (N=" << N << "; mem=" 
         << (N*(sizeof(particle_t) + sizeof(uint128_t))) << " bytes)... " << endl;
    //--------------------------------------------------------------------------
    double t0 = CPUTIME;
#define READ(_p, _n)                    \
do {                                    \
    double t = CPUTIME;                 \
    uint32_t last_i=0;                  \
    for (i=0; i < _n; i++, j++) {       \
        in >> _p;                       \
        ps[j].x[0] = _p.pos[0];         \
        ps[j].x[1] = _p.pos[1];         \
        ps[j].x[2] = _p.pos[2];         \
        ps[j].v[0] = _p.vel[0];         \
        ps[j].v[1] = _p.vel[1];         \
        ps[j].v[2] = _p.vel[2];         \
        if ((CPUTIME)-t > 2) {         \
            t = CPUTIME; \
            cerr << "\r" << ((j*100.0)/N) << "%" \
                 << " (" << ((i-last_i)/2) << "/s)                  ";\
            last_i = i; }} \
    if ((CPUTIME)-t > .02) {         \
        t = CPUTIME; cerr << "\r" << ((j*100.0)/N) << "%" << "                  \n"; } \
} while (0)

    j=0;
    READ(g, h.h_nSph);
    READ(d, h.h_nDark);
    READ(s, h.h_nStar);

    in.close();

    double t1 = CPUTIME;
    cout << "Load time: " << (t1-t0) << "s" << endl;

    //--------------------------------------------------------------------------
    cout << endl << "Pass 1 (Find range)... "; cout.flush();
    //--------------------------------------------------------------------------
    double start_time = CPUTIME;
    t0 = start_time;

    float xrange[6] = {FLT_MAX,FLT_MAX,FLT_MAX, FLT_MIN,FLT_MIN,FLT_MIN}, 
          vrange[6] = {FLT_MAX,FLT_MAX,FLT_MAX, FLT_MIN,FLT_MIN,FLT_MIN};
    for (i=0; i < N; i++)
    {
        if (ps[i].x[0] < xrange[0]) xrange[0] = ps[i].x[0];
        if (ps[i].x[1] < xrange[1]) xrange[1] = ps[i].x[1];
        if (ps[i].x[2] < xrange[2]) xrange[2] = ps[i].x[2];
        if (ps[i].x[0] > xrange[3]) xrange[3] = ps[i].x[0];
        if (ps[i].x[1] > xrange[4]) xrange[4] = ps[i].x[1];
        if (ps[i].x[2] > xrange[5]) xrange[5] = ps[i].x[2];

        if (ps[i].v[0] < vrange[0]) vrange[0] = ps[i].v[0];
        if (ps[i].v[1] < vrange[1]) vrange[1] = ps[i].v[1];
        if (ps[i].v[2] < vrange[2]) vrange[2] = ps[i].v[2];
        if (ps[i].v[0] > vrange[3]) vrange[3] = ps[i].v[0];
        if (ps[i].v[1] > vrange[4]) vrange[4] = ps[i].v[1];
        if (ps[i].v[2] > vrange[5]) vrange[5] = ps[i].v[2];
    }

    float max_xrange = max( xrange[3]-xrange[0], 
                       max( xrange[4]-xrange[1],
                            xrange[5]-xrange[2] ));

    float max_vrange = max( vrange[3]-vrange[0], 
                       max( vrange[4]-vrange[1],
                            vrange[5]-vrange[2] ));

    uint32_t nxcells = 1000;
    uint32_t nvcells = 1000;

    if (dx != 0) nxcells = (uint128_t)ceil(max_xrange / dx);
    if (dv != 0) nvcells = (uint128_t)ceil(max_xrange / dv);

    if (dx == 0) dx = max_xrange * 1.0001 / nxcells;
    if (dv == 0) dv = max_vrange * 1.0001 / nvcells;

    t1 = CPUTIME;
    cout << (t1-t0) << "s" << endl;


    //--------------------------------------------------------------------------
    cout << "Pass 2 (Generate keys; dx=" << dx 
         << " dv=" << dv 
         << " nxcells=" << nxcells
         << " nvcells=" << nvcells
         << ")... "; cout.flush();
    //--------------------------------------------------------------------------
    t0 = CPUTIME;

    uint128_t ncells1 = nxcells;
    uint128_t ncells2 = ncells1 * nxcells;
    uint128_t ncells3 = ncells2 * nxcells;
    uint128_t ncells4 = ncells3 * nvcells;
    uint128_t ncells5 = ncells4 * nvcells;

    for (i=0; i < N; i++)
    {
        keys[i] = ((uint128_t)((ps[i].x[0] - xrange[0]) / dx))
                + ((uint128_t)((ps[i].x[1] - xrange[1]) / dx)) * ncells1
                + ((uint128_t)((ps[i].x[2] - xrange[2]) / dx)) * ncells2
                + ((uint128_t)((ps[i].v[0] - vrange[0]) / dv)) * ncells3
                + ((uint128_t)((ps[i].v[1] - vrange[1]) / dv)) * ncells4
                + ((uint128_t)((ps[i].v[2] - vrange[2]) / dv)) * ncells5
                ;

        //cout << ps[i].key << endl;
        //fprintf(stdout, "%ld %f\n", (long)ps[i].key, (ps[i].x[0] - xrange[0]) / dx);
    }
    t1 = CPUTIME;
    cout << (t1-t0) << "s" << endl;
    
    //--------------------------------------------------------------------------
    cout << "Pass 3 (Sort)... "; cout.flush();
    //--------------------------------------------------------------------------
    t0 = CPUTIME;
    qsort(keys, N, sizeof(particle_t), _key_compar);
    t1 = CPUTIME;
    cout << (t1-t0) << "s" << endl;


    //--------------------------------------------------------------------------
    cout << "Pass 4 (Integrate)... "; cout.flush();
    //--------------------------------------------------------------------------
    t0 = CPUTIME;

    double    S=0;
    uint32_t  f=0;
    uint128_t key = -1;
    for (i=0; i < N; i++)
    {
        if (key != keys[i])
        {
            if (f != 0)
            {
                double f0 = ((double)f) / (double)N;
                S += f0 * log(f0);
            }

            f = 0;
            key = keys[i];
        }

        f++;
    }

    if (f != 0)
    {
        double f0 = ((double)f) / (double)N;
        S += f0 * log(f0);
    }

    t1 = CPUTIME;
    cout << (t1-t0) << "s" << endl;

    double end_time = CPUTIME;
    cout << "Compute time: " << (end_time - start_time) << "s" << endl;

    S = -S;
    cout << endl << "S=" << S << endl;

    return 0;
}

