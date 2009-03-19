#include <stdlib.h>
#include <assert.h>
#include <string>
#include <math.h>
#include <iostream>
#include "ftipsy.hpp"

using namespace std;

class Out
{
public:
    TipsyHeader h;
    ofTipsy of;
};

void help()
{
    cerr << "Usage: tipsplit -d nx,ny,nz [-std | -nat] <tipsy-in>" << endl;
    exit(2);
}

int main(int argc, char **argv)
{
    ifTipsy in;
    Out *out;
    TipsyHeader h;
    TipsyDarkParticle d;
    TipsyStarParticle s;
    TipsyGasParticle  g;
    string mode = "standard";
    string input_filename, output_prefix;

    if (argc < 2) help();

    int nx,ny,nz;
    nx = 5;
    ny = 5;
    nz = 5;

    input_filename = argv[1];
    output_prefix  = argv[2];

    in.open(input_filename.c_str(), mode.c_str());
    if (!in.is_open())
    {
        cerr << "Can't open input file " << input_filename << endl;
        exit(1);
    }

    int nCells = nx * ny * nz;
    out = new Out[nCells];

    char *output_filename = new char[input_filename.length() + 100];

    int i=0;
    for (int z=0; z < nz; z++)
        for (int y=0; y < ny; y++)
            for (int x=0; x < nx; x++, i++)
            {
                sprintf(output_filename, "%s.%i.%i.%i", output_prefix.c_str(), x,y,z);
                out[i].of.open(output_filename, mode.c_str());
                if (!out[i].of.is_open())
                {
                    cerr << "Can't open output file " << output_filename << endl;
                    exit(1);
                }
            }
    delete[] output_filename;

    in >> h;
    for (i=0; i < nCells; i++)
        out[i].of << h;
    
#define SPLIT(p, hdr_type, hdr_offs) \
    if (h.hdr_type > 0) { \
        for (i=0; i < nCells; i++) \
            out[i].of.seekp(tipsypos(tipsypos :: hdr_offs, 0)); \
        for (i=0; i < h.hdr_type; i++) { \
            in >> p; \
            int x = (int)round(nx * (p.pos[0]+0.5)); if (x == nx) x = nx-1; \
            int y = (int)round(ny * (p.pos[1]+0.5)); if (y == ny) y = ny-1; \
            int z = (int)round(nz * (p.pos[2]+0.5)); if (z == nz) z = nz-1; \
            int offs = (z*nx*ny) + (y*nx) + x; \
            out[offs].of << p; \
            out[offs].h.hdr_type++; }}


    SPLIT(g, h_nGas,  gas);
    SPLIT(d, h_nDark, dark);
    SPLIT(s, h_nStar, star);

    in.close();
    for (i=0; i < nCells; i++)
    {
        out[i].of.seekp(tipsypos(tipsypos::header, 0));
        out[i].h.h_time  = h.h_time;
        out[i].h.h_nDims = h.h_nDims;
        out[i].h.h_nBodies = out[i].h.h_nGas + out[i].h.h_nDark + out[i].h.h_nStar;
        out[i].of << out[i].h;
        out[i].of.close();
    }

    return 0;
}

