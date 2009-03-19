#include <stdlib.h>
#include <assert.h>
#include <string>
#include <math.h>
#include <complex>
#include <iostream>
#include <pthread.h>
#include "ftipsy.hpp"

using namespace std;

void help()
{
    cerr << "Usage: cutbox [--no-wrap] [-std | -nat] <tipsy-in> <tipsy-out> x y z rx ry rz" << endl;
    exit(2);
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
    bool force_no_wrap = false;
    bool center = false;

    if (argc < 7) help();

    argv++; // skip argv[0]
    for (i=0; i < argc; i++, argv++)
    {
        if (!strcmp("--no-wrap", *argv))
            force_no_wrap = true;
        else if (!strcmp("-std", *argv))
            {} /* already default */
        else if (!strcmp("-nat", *argv))
            mode = "native";
        else if (!strcmp("--center", *argv))
            center = true;
        else
            break;
    }

    if (argc-i < 7) help();

    string input_filename(*argv++);
    string output_filename(*argv++);

    double cx = atof(*argv++);
    double cy = atof(*argv++);
    double cz = atof(*argv++);
    double rx = atof(*argv++);

    double ry = -1; 
    double rz = -1;

    if (argc-i > 7)
    {
        ry = atof(*argv++);
        rz = atof(*argv++);
    }

    double left   = cx-rx;
    double right  = cx+rx;
    double top    = cy-ry;
    double bottom = cy+ry;
    double back   = cz-rz;
    double front  = cz+rz;

    const bool wraplr = (left < -0.5) || (right  > 0.5);
    const bool wraptb = (top  < -0.5) || (bottom > 0.5);
    const bool wrapbf = (back < -0.5) || (front  > 0.5);

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

    oh.h_nBodies = oh.h_nGas = oh.h_nDark = oh.h_nStar = 0;

    if (force_no_wrap || (!wraplr && !wraptb && !wrapbf))
    {
        for (i=0; i < h.h_nDark; i++)
        {
            in >> d;
            const double x = d.pos[0];
            const double y = d.pos[1];
            const double z = d.pos[2];

            bool accept;
            
            if (ry != -1)
            {
                accept =  (left <= x&&x <= right)
                       && (top  <= y&&y <= bottom)
                       && (back <= z&&z <= front);
            }
            else
            {
                accept = sqrt(pow(x-cx,2) + pow(y-cy,2) + pow(z-cz,2)) <= rx;
            }

            if (accept) 
            {
                if (center)
                {
                    d.pos[0] = x - cx;
                    d.pos[1] = y - cy;
                    d.pos[2] = z - cz;
                }
                out << d;
                oh.h_nDark++;
            }
        }
    }
    else
    {
        cerr << "Wrapping." << endl;
        cerr << "BOX BEFORE: ";
        cerr << "left "   << left   << "  ";
        cerr << "right "  << right  << "  ";
        cerr << "top "    << top    << "  ";
        cerr << "bottom " << bottom << "  ";
        cerr << "front "  << front  << "  ";
        cerr << "back "   << back   << endl;

        left   = fmod(left   + 1.5, 1.0) - 0.5;
        right  = fmod(right  + 1.5, 1.0) - 0.5;
        top    = fmod(top    + 1.5, 1.0) - 0.5;
        bottom = fmod(bottom + 1.5, 1.0) - 0.5;
        front  = fmod(front  + 1.5, 1.0) - 0.5;
        back   = fmod(back   + 1.5, 1.0) - 0.5;

        cerr << "BOX AFTER:  ";
        cerr << "left "   << left   << "  ";
        cerr << "right "  << right  << "  ";
        cerr << "top "    << top    << "  ";
        cerr << "bottom " << bottom << "  ";
        cerr << "front "  << front  << "  ";
        cerr << "back "   << back   << endl;

        for (i=0; i < h.h_nDark; i++)
        {
            in >> d;
            const double x = d.pos[0];
            const double y = d.pos[1];
            const double z = d.pos[2];
            bool accept =  (wraplr ? (left <= x || x <= right)  : (left <= x&&x <= right))
                        && (wraptb ? (top  <= y || y <= bottom) : (top  <= y&&y <= bottom))
                        && (wrapbf ? (back <= z || z <= front)  : (back <= z&&z <= front));

            if (accept) 
            {
                d.pos[0] = fmod(x - cx + 1.5, 1.0) - 0.5;
                d.pos[1] = fmod(y - cy + 1.5, 1.0) - 0.5;
                d.pos[2] = fmod(z - cz + 1.5, 1.0) - 0.5;
                out << d;
                oh.h_nDark++;
            }
        }
    }

    out.seekp(tipsypos(tipsypos::header, 0));

    oh.h_nDims = h.h_nDims;
    oh.h_time = h.h_time;
    oh.h_nBodies = oh.h_nGas + oh.h_nDark + oh.h_nStar;
    out << oh;

    in.close();
    out.close();

    return 0;
}

