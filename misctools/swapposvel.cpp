#include <stdlib.h>
#include <assert.h>
#include <string>
#include <iostream>
#include <math.h>
#include "ftipsy.hpp"

using namespace std;

void help()
{
    cerr << "Usage: swapposvel [-std | -nat] <tipsy-in> <tipsy-out>" << endl;
    exit(2);
}

void doswap(TipsyBaseParticle &p)
{
    const float x = p.pos[0];
    const float y = p.pos[1];
    const float z = p.pos[2];

    p.pos[0] = p.vel[0] * 22.0;
    p.pos[1] = p.vel[1] * 22.0;
    p.pos[2] = p.vel[2] * 22.0;

    p.vel[0] = x;
    p.vel[1] = y;
    p.vel[2] = z;
}

int main(int argc, char **argv)
{
    int i;
    ifTipsy in;
    ofTipsy out;
    TipsyHeader h;
    TipsyHeader oh;
    TipsyDarkParticle d;
    TipsyStarParticle s;
    TipsyGasParticle  g;
    string mode = "standard";

    if (argc < 3) help();

    argv++; // skip argv[0]
    for (i=0; i < argc; i++, argv++)
    {
        if (!strcmp("-std", *argv))
            {} /* already default */
        else if (!strcmp("-nat", *argv))
            mode = "native";
        else
            break;
    }

    if (argc-i < 3) help();

    string input_filename(*argv++);
    string output_filename(*argv++);

    in.open(input_filename.c_str(), mode.c_str());
    if (!in.is_open())
    {
        cerr << "Can't open input file " << input_filename << endl;
        exit(1);
    }

    out.open(output_filename.c_str(), mode.c_str());
    if (!out.is_open())
    {
        cerr << "Can't open output file " << output_filename << endl;
        exit(1);
    }

    in  >> h;
    out << h;

#define SWAP(p, N) \
    for (i=0; i < N; i++) { in >> p; doswap(p); out << p; }

    SWAP(g, h.h_nSph);
    SWAP(d, h.h_nDark);
    SWAP(s, h.h_nStar);

    in.close();
    out.close();

    return 0;
}

