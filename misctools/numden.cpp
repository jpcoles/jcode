#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "ftipsy.hpp"

using namespace std;

void help()
{
    cerr << "Usage: numden tipsy xdim ydim zdim" << endl;
    exit(2);
}

int main(int argc, char **argv)
{
    ifTipsy in;
    TipsyHeader h;
    TipsyDarkParticle d;

    if (argc < 4) help();

    in.open(argv[1], "standard");
    if (!in.is_open())
    {
        cerr << "Unable to open input file " << argv[1] << endl;
        exit(2);
    }

    const int xdim = atoi(argv[2]);
    const int ydim = atoi(argv[3]);
    const int zdim = atoi(argv[4]);

    h.h_time = 0;
    h.h_nDims = 3;
    h.h_nBodies = h.h_nDark = xdim * ydim * zdim;
    h.h_nStar = h.h_nGas = 0;

    in >> h;

    uint32_t nCells = xdim * ydim * zdim;
    uint32_t *den = new uint32_t[nCells];
    memset(den, 0, nCells * sizeof(uint32_t));

    for (uint32_t i=0; i < h.h_nBodies; i++)
    {
        in >> d;
        int x = (d.pos[0] + 0.5) * xdim;
        int y = (d.pos[1] + 0.5) * ydim;
        int z = (d.pos[2] + 0.5) * zdim;
        den[z * (ydim*xdim) + y*xdim + x]++;
    }

    in.close();

    printf("%i\n", nCells);
    for (uint32_t i=0; i < nCells; i++)
        printf("%24.15e\n", (double)den[i] / nCells);

    return 0;
}

