#include <stdlib.h>
#include <assert.h>
#include <string>
#include <iostream>
#include "ftipsy.hpp"

using namespace std;

void help()
{
    cerr << "Usage: tipexplode [-std | -nat] <tipsy-in> <tipsy-out> <factor>" << endl;
    exit(2);
}

void explode(TipsyBaseParticle &p, float fac)
{
    p.pos[0] += p.vel[0] * fac;
    p.pos[1] += p.vel[1] * fac;
    p.pos[2] += p.vel[2] * fac;
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

    if (argc < 4) help();

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

    if (argc-i < 4) help();

    string input_filename(*argv++);
    string output_filename(*argv++);
    float  fac = atof(*argv++);

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

#define EXPLODE(p, N) \
    for (i=0; i < N; i++) { in >> p; explode(p, fac); out << p; }

    EXPLODE(g, h.h_nSph);
    EXPLODE(d, h.h_nDark);
    EXPLODE(s, h.h_nStar);

    in.close();
    out.close();

    return 0;
}

