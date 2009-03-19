#include <stdlib.h>
#include <assert.h>
#include <string>
#include <math.h>
#include <iostream>
#include "ftipsy.hpp"

using namespace std;

// Take every nth particle of the given type.
const int gasFrac  = 10000;
const int darkFrac = 10000;
const int starFrac = 10000;

void help()
{
    cerr << "Usage: tipsubsampl <tipsy-in> <tipsy-out>" << endl;
    cerr << "Produces a new tipsy file every nth particle is taken from the original file." << endl;
    cerr << "n_gas=" << gasFrac << " ";
    cerr << "n_dark=" << darkFrac << " ";
    cerr << "n_star=" << starFrac << endl;
    exit(2);
}

int main(int argc, char **argv)
{
    ifTipsy in;
    ofTipsy out;
    TipsyHeader ih, oh;
    TipsyDarkParticle d;
    TipsyStarParticle s;
    TipsyGasParticle  g;
    string mode = "standard";
    string input_filename, output_filename;

    if (argc < 3) help();

    input_filename  = argv[1];
    output_filename = argv[2];

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

    in >> ih;
    memcpy(&oh, &ih, sizeof(oh));


    out << oh;

    oh.h_nBodies = oh.h_nGas = oh.h_nDark = oh.h_nStar = 0;

    cerr << "gasFrac  = " << gasFrac  << endl;
    cerr << "darkFrac = " << darkFrac << endl;
    cerr << "starFrac = " << starFrac << endl;

#define SAMPLE1(p, hdr_type, hdr_offs, frac) \
    for (int i=0; i < ih.hdr_type; i+=frac) { \
        in.seekg(tipsypos(tipsypos::hdr_offs, i)); \
        in >> p; out << p; \
        oh.hdr_type++; }

#define SAMPLE2(p, hdr_type, hdr_offs, frac) \
    for (int i=0; i < ih.hdr_type; i+=frac) { \
        long r = random() % frac; \
        if (i+r >= ih.hdr_type) break; \
        in.seekg(tipsypos(tipsypos::hdr_offs, i+r)); \
        in >> p; out << p; \
        oh.hdr_type++; }

#define SAMPLE SAMPLE2

    SAMPLE(g, h_nGas,  gas,  gasFrac);
    SAMPLE(d, h_nDark, dark, darkFrac);
    SAMPLE(s, h_nStar, star, starFrac);

    oh.h_nBodies = oh.h_nGas + oh.h_nDark + oh.h_nStar;
    cerr << "nBodies = " << oh.h_nBodies << endl;

    out.seekp(tipsypos(tipsypos::header, 0));
    out << oh;

    in.close();
    out.close();

    return 0;
}

