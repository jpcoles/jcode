#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "ftipsy.hpp"

using namespace std;

void help()
{
    cerr << "Usage: minpot tipsy" << endl;
    exit(2);
}

int main(int argc, char **argv)
{
    ifTipsy in;
    TipsyHeader h;
    TipsyDarkParticle d;
    TipsyDarkParticle dmin, dmax;

    if (argc < 2) help();

    in.open(argv[1], "standard");
    if (!in.is_open())
    {
        cerr << "Unable to open input file " << argv[1] << endl;
        exit(2);
    }

    in >> h;

    dmin.phi = 1e30;
    dmax.phi = -1e30;

    uint32_t min_i=0;
    uint32_t max_i=0;
    for (uint32_t i=0; i < h.h_nBodies; i++)
    {
        in >> d;
        if (d.phi < dmin.phi) { dmin = d; min_i = i+1; }
        if (d.phi > dmax.phi) { dmax = d; max_i = i+1; }
    }

    in.close();

    if (min_i == 0 || max_i == 0)
    {
        fprintf(stderr, "Couldn't find min/max\n");
        exit(1);
    }
    else
    {
        printf("%.4f ", 1.0/h.h_time - 1);
        printf("%i %.7e %.7e %.7e %.7e ",  min_i, (dmin.pos[0]+0.5)*65.7, (dmin.pos[1]+0.5)*65.7, (dmin.pos[2]+0.5)*65.7, dmin.phi);
        printf("%i %.7e %.7e %.7e %.7e\n", max_i, (dmax.pos[0]+0.5)*65.7, (dmax.pos[1]+0.5)*65.7, (dmax.pos[2]+0.5)*65.7, dmax.phi);
    }

    return 0;
}

