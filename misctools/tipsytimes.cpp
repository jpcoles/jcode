#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <assert.h>
#include "ftipsy.hpp"

using namespace std;

int precision=4;

typedef void (*output_func)(bool, bool, int, int, float);

void normal_output(bool output_time, bool output_z, int mod, int nfiles, float time)
{
    if (output_time) printf("%10.*f ", precision, time);
    if (output_z)    printf("%10.*f ", precision, 1/time - 1);
    printf("\n");
}

void python_output(bool output_time, bool output_z, int mod, int nfiles, float time)
{
    if (nfiles != 1) printf(", ");
    if (nfiles % mod == 0) printf("\n");
    if (output_time) printf("'%0.*f'", precision, time);
    if (output_z)    printf("'%0.*f'", precision, 1/time - 1);
}

void bash_output(bool output_time, bool output_z, int mod, int nfiles, float time)
{
    if (nfiles != 1) printf(" ");
    if (nfiles % mod == 0) printf("\n    ");
    if (output_time) printf("%0.*f", precision, time);
    if (output_z)    printf("%0.*f", precision, 1/time - 1);
}

void output(int argc, char **argv, int argind, int output_time, int output_z, output_func f, int mod)
{
    for (int nfiles=0; argind < argc; argind++)
    {
        ifTipsy in;
        TipsyHeader h;

        in.open(argv[argind], "standard");
        if (!in.is_open())
        {
            fprintf(stderr, "ERROR: Can't open %s\n", argv[argind]);
            continue;
        }

        nfiles++;

        in >> h;

        f(output_time, output_z, mod, nfiles, h.h_time);

        in.close();
    }
}

void help()
{
    fprintf(stderr, 
        "Usage: tipsytimes [-t] [-z] [-d #] [--python] tipsy-files\n"
        "\n"
        "    -t                     Display the time\n"
        "    -z                     Display the redshift\n"
        "    --python               Output python lists\n"
        "    --bash                 Output bash lists\n"
        "    -d #                   Display fractions with # digits\n"
        "\n"
        "If --python is not specified and both -t -z are, then two\n"
        "columns are displayed with the first column being the time\n"
        "and the second being the redshift\n"
        "\n");
    exit(2);
}

int main(int argc, char **argv)
{
    bool output_time = false;
    bool output_z = false;
    bool output_python = false;
    bool output_bash = false;

    int mod=0;

    while (1) 
    {
        int c;
        //int this_option_optind = optind ? optind : 1;
        int option_index = 0;
        static struct option long_options[] = {
            {"python", 0, 0, 0},
            {"bash", 0, 0, 0},
            {"help", 0, 0, 'h'},
            {0, 0, 0, 0}
        };

        c = getopt_long(argc, argv, "tzd:h", long_options, &option_index);
        if (c == -1)
            break;

        switch (c)
        {
            case 0:
                if (!strcmp("python", long_options[option_index].name))
                    output_python = true;
                else if (!strcmp("bash", long_options[option_index].name))
                    output_bash = true;
                break;
            case 't':
                output_time = true;
                break;
            case 'z':
                output_z = true;
                break;
            case 'd':
                precision = atoi(optarg);
                if (precision == 0)
                {
                    fprintf(stderr, "ERROR: Bad number of digits %s\n", optarg);
                    exit(2);
                }
                break;

            case 'h':
                help();

            default:
                break;
        }
    }

    if (! (output_time || output_z) )
        output_time = true;

    if (optind == argc)
        help();

    if (output_python)
    {
        mod = 80 / (precision+5);
        if (mod < 1) mod = 1;

        if (output_time)
        {
            printf("tipsy_times = [");
            output(argc, argv, optind, output_time, false, python_output, mod);
            printf("]\n\n");
        }
        if (output_z)
        {
            printf("tipsy_z = [");
            output(argc, argv, optind, false, output_z, python_output, mod);
            printf("]\n\n");
        }
    }
    else if (output_bash)
    {
        mod = 80 / (precision+5);
        if (mod < 1) mod = 1;

        if (output_time)
        {
            printf("TIPSY_TIMES = '");
            output(argc, argv, optind, output_time, false, bash_output, mod);
            printf("'\n\n");
        }
        if (output_z)
        {
            printf("TIPSY_Z = '");
            output(argc, argv, optind, false, output_z, bash_output, mod);
            printf("'\n\n");
        }
    }
    else
    {
        output(argc, argv, optind, output_time, output_z, normal_output, mod);
    }
        

    return 0;

}

