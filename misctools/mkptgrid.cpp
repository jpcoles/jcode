#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "ftipsy.hpp"

using namespace std;

void help()
{
    cerr << "Usage: mkptgrid tipsyout xdim ydim zdim" << endl;
    exit(2);
}

int main(int argc, char **argv)
{
    ofTipsy out;
    TipsyHeader h;
    TipsyDarkParticle d;

    if (argc < 5) help();

    out.open(argv[1], "standard");
    if (!out.is_open())
    {
        cerr << "Unable to create output file " << argv[1] << endl;
        exit(2);
    }

    const int xdim = atoi(argv[2]);
    const int ydim = atoi(argv[3]);
    const int zdim = atoi(argv[4]);

    h.h_time = 0;
    h.h_nDims = 3;
    h.h_nBodies = h.h_nDark = xdim * ydim * zdim;
    h.h_nStar = h.h_nGas = 0;

    out << h;

    out.seekp(tipsypos(tipsypos::dark, 0));
    for (int z=0; z < zdim; z++)
        for (int y=0; y < ydim; y++)
            for (int x=0; x < xdim; x++)
            {
                d.pos[0] = x / (double)xdim - 0.5;
                d.pos[1] = y / (double)ydim - 0.5;
                d.pos[2] = z / (double)zdim - 0.5;
                out << d;
            }

    out.close();

    return 0;
}

