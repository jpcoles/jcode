/* 
 * Particle INTegrator
 *
 * Jonathan Coles and Joachim Stadel
 *
 * $Revision$
 */

#define PINT_MAIN

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <math.h>
#include <getopt.h>
#include <assert.h>

#include "io.h"
#include "pint.h"
#include "tipsy.h"

#if USE_TREE_JS
#include "tree_js.h"
#endif

#if USE_TREE_JPC
#include "tree_jpc.h"
#endif

#define VERSION "0.1"
#define VERSION_DATE "20060419"
#define CVS_REVISION "$Revision$"

real MMAX(real a, real b)
{
    if (a > b) return a;
    return b;
}

real MMIN(real a, real b)
{
    if (a < b) return a;
    return b;
}

real SQDIST(struct particle *p, struct particle *q)
{
    int i;
    real dist=0;
    for (i=0; i < 3; i++)
    {
        real d = POSAT(q, i, p->timeNext) - p->r[i];
        dist += d*d;
    }
    assert(dist >= 0);

    return dist;
}

real DIST(struct particle *p, struct particle *q)
{
    return sqrt(SQDIST(p, q));
}

real SQMAG(real *a)
{
    int i;
    real dist=0;
    for (i=0; i < 3; i++)
    {
        //double d = a[i];
        dist += a[i]*a[i];
    }

    assert(dist >= 0);

    return dist;
}

real MAG(real *a)
{
    return sqrt(SQMAG(a));
}


static void help()
{
    fprintf(err, "pint [-f input-file] [-o output-file]\n");
    exit(2);
}

int main(int argc, char **argv)
{
    int i;

    struct env env;

    in = stdin;
    out = stdout;
    err = stderr;
    logfp = NULL;

    static struct option long_options[] = {
        {"file", required_argument, 0, 'f'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    infile = NULL;
    outfilebase = NULL;
    verbosity = 0;
    memset(&env, 0, sizeof(env));

    /*========================================================================
     * Process the command line flags
     *======================================================================*/
    while (1)
    {
        int option_index = 0;
        int c = getopt_long(argc, argv, "f:vlo:",
                            long_options, &option_index);

        if (c == -1) break;

        switch (c)
        {
            case 0:
                break;

            case 'f': infile      = optarg;       break;
            case 'o': outfilebase = optarg;       break;
            case 'v': verbosity++;                break;
            case 'l': loglevel++;                 break;

            case 'h': help(); break;
            case '?': break;
        }
    }

    if (loglevel != 0)
    {
        int i;
        char logfile[256];

        snprintf(logfile, 256, "%s.log.%i", "pint", getpid());
        for (i=1; i <= 1000 && (!access(logfile, R_OK) || !access(logfile, W_OK)); i++)
        {
            snprintf(logfile, 256, "%s.log.%i-%i", "pint", getpid(), i);
        }

        logfp = fopen(logfile, "w");
        if (logfp == NULL) 
        {
            fprintf(err, "WARNING: Can't create log file %s. No log will be kept.", logfile);
        }

        VL(1) fprintf(out, "Logfile: %s\n", logfile);
    }

    VL(1) fprintf(out, "Verbosity level %i\n", verbosity);

    sleep(5);

    /*========================================================================
     * Load the initial conditions
     *======================================================================*/
    if (infile != NULL)
    {
        LOG(1) fprintf(logfp, "Loading %s.\n", infile);
        if (judge_and_load_file(infile, &env))
        {
            fprintf(err, "Unable to load file.\n");
            exit(1);
        }
    }
    else
    {
        ic_threebody(&env);
    }

    if (outfilebase != NULL)
    {
    }


    for (i=0; i < 15; i++)
    {

#if USE_TREE_JS
        tree_build_js(&env);
        tree_free_js(&env);
#endif

#if USE_TREE_JPC
        tree_build_jpc(&env);
        tree_free_jpc(&env);
#endif
    }

    free(env.ps);
    env.ps = NULL;
    free(env.p);
    env.p = NULL;

    if (logfp != NULL && logfp != stdin && logfp != in && logfp != out && logfp != err)
        fclose(logfp);

    return 0;
}

