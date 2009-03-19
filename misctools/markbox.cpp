#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string>
#include <math.h>
#include <complex>
#include <pthread.h>
#include "ftipsy.hpp"

using namespace std;

void help()
{
    fprintf(stderr, "Usage: markbox [-grp] [-std | -nat] <tipsy-in> <out> x y z r");
    exit(2);
}

int main(int argc, char **argv)
{
    int i;
    ifTipsy in;
    FILE *out;
    TipsyHeader h;
    TipsyDarkParticle d;
    TipsyStarParticle s;
    TipsyGasParticle  g;
    string mode = "standard";
    bool center = false;
    bool group  = false;

    if (argc < 7) help();

    argv++; // skip argv[0]
    for (i=0; i < argc; i++, argv++)
    {
        if (!strcmp("-grp", *argv))
            group = true;
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
    double r  = atof(*argv++);

    in.open(input_filename.c_str(), mode.c_str());
    if (!in.is_open())
    {
        fprintf(stderr, "Can't open input file %s\n", input_filename.c_str());
        exit(1);
    }

    out = fopen(output_filename.c_str(), "w");
    if (out == NULL)
    {
        fprintf(stderr, "Can't open output file %s\n", output_filename.c_str());
        exit(1);
    }

    in >> h;

    if (!group)
        fprintf(out, "%i 0 0\n", h.h_nDark);
    else
        fprintf(out, "%i\n", h.h_nDark);

    for (i=0; i < h.h_nDark; i++)
    {
        in >> d;
        const double x = d.pos[0];
        const double y = d.pos[1];
        const double z = d.pos[2];

        bool accept = sqrt(pow(x-cx,2) + pow(y-cy,2) + pow(z-cz,2)) <= r;

        if (group)
            fprintf(out, "%i\n", (int)accept);
        else
            if (accept) 
                fprintf(out, "%i\n", i+1);
    }

    in.close();
    fclose(out);

    return 0;
}

