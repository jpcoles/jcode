#include <iostream>
#include <sstream>
#include <string>
#include <math.h>
#include <assert.h>

#include "ftipsy.hpp"

using namespace std;

#define VERSION "0.1"

#define SCALE(p) \
do { \
    (p).pos[0] += dx; (p).pos[1] += dy; (p).pos[2] += dz; \
    if ((p).pos[0] > 0.5F) (p).pos[0] -= 1.0F; else if ((p).pos[0] < -0.5F) (p).pos[0] += 1.0F;\
    if ((p).pos[1] > 0.5F) (p).pos[1] -= 1.0F; else if ((p).pos[1] < -0.5F) (p).pos[1] += 1.0F;\
    if ((p).pos[2] > 0.5F) (p).pos[2] -= 1.0F; else if ((p).pos[2] < -0.5F) (p).pos[2] += 1.0F;\
    if ((pi % two_percent) == 0) cerr << "."; \
} while(0)

void help() 
{
    cerr << "Usage: tipscale <in tipsy file> <out tipsy file> dx dy dz" << endl
         << endl
         << "Send questions, comments, bug reports to Jonathan Coles <jonathan@physik.uzh.ch>" << endl
         << endl;
    
    exit(2);
}

void version()
{
    cerr << "tipscale v" << VERSION << endl
         << "Tipsy Scaler." << endl
         << "Written by Jonathan Coles <jonathan@physik.uzh.ch>" << endl;

    exit(0);
}

int main(int argc, char **argv)
{
    int i;
    bool is_std = true;
    ifTipsy in;
    ofTipsy out;
    TipsyHeader h;
    TipsyGasParticle  g; // A gas particle
    TipsyDarkParticle d; // A dark particle
    TipsyStarParticle s; // A star particle

    int optind = 1;
    if (argc >= 2 && !strcmp(argv[1], "--help")) help();
    if (argc >= 2 && !strcmp(argv[1], "--version")) version();
    if (argc < 6) help();

    in.open(argv[optind], is_std ? "standard" : "native");
    if (!in.is_open()) 
    {
        cerr << "ERROR: Unable to open Tipsy binary " << argv[optind] << endl;
        exit(1);
    }
    optind++;

    out.open(argv[optind], is_std ? "standard" : "native");
    if (!in.is_open()) 
    {
        cerr << "ERROR: Unable to open Tipsy binary " << argv[optind] << endl;
        exit(1);
    }
    optind++;

    float dx = atof(argv[optind++]);
    float dy = atof(argv[optind++]);
    float dz = atof(argv[optind++]);

    in >> h; out << h;

    cerr << "   ""                                                   ] 100%\r";
    cerr << "0% [";

    unsigned int two_percent = (int)(h.h_nBodies * 0.02);
    unsigned long pi=0;
    for( i=0; i<h.h_nSph;  i++, pi++ ) { in >> g; SCALE(g); out << g; }
    for( i=0; i<h.h_nDark; i++, pi++ ) { in >> d; SCALE(d); out << d; }
    for( i=0; i<h.h_nStar; i++, pi++ ) { in >> s; SCALE(s); out << s; }


    in.close();
    out.close();


    return 0;
}

