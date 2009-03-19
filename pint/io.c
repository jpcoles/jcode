#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include <assert.h>
#include <time.h>
#include "pint.h"
#include "io.h"

#ifndef M_PI
#define M_PI 3.1415926535
#endif

extern int AccelCount;

void print_queue_helper(struct particle *p, int lvl)
{
    if (p == NULL) return;
    fprintf(err, "\n%*sid: %i  d(%i) r(%e %e %e) v(%e %e %e) a(%e %e %e)\n", // t(%e %e)\n", 
        lvl*4,
        "", 
        ID(p), 
        p->pqDist,
        p->r[0], p->r[1], p->r[2],
        p->v[0], p->v[1], p->v[2],
        p->a[0], p->a[1], p->a[2]
        //p->time, p->timeNext
        );
    //fprintf(err, "L: "); 
    //print_queue_helper(p->pqLeft, lvl+1);
    //fprintf(err, "R: "); 
    //print_queue_helper(p->pqRight, lvl+1);
}

void print_queue(struct particle *p)
{
    print_queue_helper(p, 0);
}

/*============================================================================
 *                               loadascii
 *
 * Load initial conditions from an ascii file in the format
 *      nparticles
 *      ptype pmass x y z vx vy vz
 *      ptype pmass x y z vx vy vz
 *         ... 
 *==========================================================================*/
int loadascii(char *filename, struct env *env)
{
    int i;

    if (verbosity > 0)
        fprintf(err, "Loading ascii format.\n");

    FILE *fp = fopen(filename, "r");
    if (fp == NULL) return 1;

    fscanf(fp, "%i\n", &(env->N));
    double rx, ry, rz,
          vx, vy, vz,
          m;

    int t;

    env->p = (struct particle **)malloc(env->N * sizeof(struct particle *));
    if (env->p == NULL) return 1;

    struct particle *tp = (struct particle *)calloc(env->N, sizeof(struct particle));
    if (tp == NULL) return 1;

    for (i=0; i < env->N; i++)
    {
        fscanf(fp, "%i %lf %lf %lf %lf %lf %lf %lf\n", &t, &m, &rx, &ry, &rz, &vx, &vy, &vz);

        env->p[i] = &(tp[i]); //(struct particle *)calloc(1, sizeof(struct particle));
        PQ_INIT(env->p[i]);
        ID(env->p[i]) = i;
        SOFT(env->p[i]) = DEFAULT_SOFTENING;
        MASS(env->p[i]) = m;
        env->p[i]->r[0] = rx;
        env->p[i]->r[1] = ry;
        env->p[i]->r[2] = rz;
        env->p[i]->v[0] = vx;
        env->p[i]->v[1] = vy;
        env->p[i]->v[2] = vz;

        fprintf(err, "%i %f %f %f %f %f %f %f\n", t, m, rx, ry, rz, vx, vy, vz);
    }

    fclose(fp);

    return 0;
}

/*============================================================================
 *                               loadtipsy
 *
 * Use the given Tipsy context to load particles.
 *==========================================================================*/
int loadtipsy(TCTX tctx, struct env *env)
{
    int i, ret;

    if (verbosity > 0)
        fprintf(err, "Loading tipsy format (%s).\n",
            tctx->bNative ? "native" : "non-native");

    if ((ret = TipsyReadAll(tctx))) return ret;

    env->N = tctx->nGas + tctx->nDark + tctx->nStar;
    env->p = (struct particle **)malloc(env->N * sizeof(struct particle *));
    if (env->p == NULL) return 1;

    env->ps = (struct particle *)calloc(env->N, sizeof(struct particle));
    if (env->ps == NULL) return 1;

    for (i=0; i < env->N; i++)
    {
        int j;
        int type;
        double soft;
        struct base_particle *p = pTipsyParticle(tctx, i, &type, &soft);
        assert(p != NULL);

        env->p[i] = &(env->ps[i]); //(struct particle *)calloc(1, sizeof(struct particle));

        PQ_INIT(env->p[i]);
        ID(env->p[i]) = i;
        SOFT(env->p[i]) = DEFAULT_SOFTENING;

        for (j=0; j < 3; j++) env->p[i]->r[j] = p->pos[j];
        for (j=0; j < 3; j++) env->p[i]->v[j] = p->vel[j];
        MASS(env->p[i]) = p->mass;
    }

    TipsyFinish(tctx);
    return 0;
}

/*============================================================================
 *                          judge_and_load_file
 *
 * Attempt to determine the file type and if successful, load the file.
 * This routine will recognize tipsy format (native and non-native) and an
 * ascii format (see loadascii).
 *
 * Returns 0 if the file was successfully loaded and 1 otherwise.
 *==========================================================================*/
int judge_and_load_file(char *filename, struct env *env)
{
    TCTX tctx = NULL;

    /* 1) try tipsy non-native format */
    if (TipsyInitialize(&tctx, 0, filename) || loadtipsy(tctx, env))
        /* 2) try tipsy native format */
        if (TipsyInitialize(&tctx, 1, filename) || loadtipsy(tctx, env))
            /* 3) try ascii file */
            if (loadascii(filename, env))
            {
                fprintf(err, "File %s is in an unsupported format.\n", filename);
                return 1;
            }

    return 0;
}

/*============================================================================
 *                                write_tipsy
 *
 * Write the current particle data to a file in tipsy non-native format.
 * The filename used is in the form <filename>-<seqno>.dat, where <seqno> is
 * zero-padded to five places.
 *==========================================================================*/
void write_tipsy(char *filename, int seqno, struct env *env)
{
    int i, j;
    char name[20];
    snprintf(name, 19, "%s-%05i.dat", filename, seqno);

    TCTX tctx;
    TipsyInitialize(&tctx, 0, NULL);

    for (i=0; i < env->N; i++)
    {
        struct star_particle *p = (struct star_particle *)calloc(1, sizeof(struct star_particle));
        p->mass = MASS(ps[i]);
        for (j=0; j < 3; j++) p->pos[j] = ps[i]->r[j];
        for (j=0; j < 3; j++) p->vel[j] = ps[i]->v[j];
        TipsyAddStar(tctx, p);
    }

    TipsyWriteAll(tctx, 0, name);

    TipsyFinish(tctx);
}

void write_ascii_fp(FILE *fp, struct env *env)
{
    int i;
    for (i=0; i < env->N; i++)
    {
        fprintf(fp, "%i %e %e %e %e %e %e %e\n", 
            0, MASS(env->p[i]), 
            env->p[i]->r[0], env->p[i]->r[1], env->p[i]->r[2], 
            env->p[i]->v[0], env->p[i]->v[1], env->p[i]->v[2]);
    }
    fprintf(fp, "\n\n");
}

void write_ascii(char *filename, int seqno, struct env *env)
{
    char name[23];
    snprintf(name, 22, "%s-%05i.asc.dat", filename, seqno);

    FILE *fp = fopen(name, "w");
    if (fp == NULL) return;

    fprintf(fp, "%i\n", env->N);
    write_ascii_fp(fp, env);
    fclose(fp);
}

void write_density(FILE *fp, struct env *env)
{
    int i;
    double bins[100];

    double rmax = MAG(POS(env->p[env->N-1]));

    double db = pow(rmax, 0.01);
    double binmin = 0;
    double binmax = db;
    int b = 0;

    for (i=0; i < 100; i++) bins[i] = 0;
    
    for (i=0; i < env->N; i++)
    {
        if (MAG(POS(env->p[i])) > binmax) 
        {
            bins[b] /= 4.*M_PI*pow(binmax, 2) * (binmax-binmin);
            if (b == 0) bins[0] /= 3.;

            fprintf(fp, "%22.15e %22.15e\n", binmax, bins[b]);
            binmin = binmax;
            binmax *= db;
            b++;
        }

        bins[b] += MASS(env->p[i]);
    }

    fprintf(fp, "\n\n");
    fflush(fp);
    //fprintf(err, "density profile written.\n");
}

void write_timesteps(FILE *fp, struct env *env)
{
#if 0
    int i;
    double part = .10;
    double mass = 0;
    fprintf(fp, "%22.15e", time);
    for (i=0; i < env->N && part <= 1; i++)
    {
        mass += MASS(env->p[i]);
        while (mass > part*TOTAL_MASS && part <= 1)
        {
            fprintf(fp, " %e", MAG(RADIUS(env->p[i])));
            part += .10;
        }
    }
    fprintf(fp, "\n");
    fflush(fp);
#endif
}

#if 0
void write_energy(FILE *fp, double time)
{
    double E, T, V;
    E = TOTAL_ENERGY(env->p, &T, &V);
    //fprintf(fp, "%22.15e %22.15e %22.15e %22.15e\n", time, pow(E*5e3, 2), T, V);
    fprintf(fp, "%22.15e %22.15e %22.15e %22.15e\n", time, E, T, V);
    fflush(fp);
}

void write_mass(FILE *fp, double time)
{
    int i;
    double part = .10;
    double mass = 0;
    fprintf(fp, "%22.15e", time);
    for (i=0; i < env->N && part <= 1; i++)
    {
        mass += MASS(env->p[i]);
        while (mass > part*TOTAL_MASS && part <= 1)
        {
            fprintf(fp, " %e", MAG(RADIUS(env->p[i])));
            part += .10;
        }
    }
    fprintf(fp, "\n");
    fflush(fp);
}
#endif

void ic_binarystar(struct env *env)
{
    int i, j;

    if (verbosity > 0) fprintf(err, "Using binary star system.\n");

    env->N = 2;
    env->p = (struct particle **)malloc(env->N * sizeof(struct particle *));
    double cm[3];

    env->p[0] = (struct particle *)calloc(1, sizeof(struct particle));
    PQ_INIT(env->p[0]);
    ID(env->p[0]) = 0;
    MASS(env->p[0]) = 0.5;
    SOFT(env->p[0]) = DEFAULT_SOFTENING;

    env->p[1] = (struct particle *)calloc(1, sizeof(struct particle));
    PQ_INIT(env->p[1]);
    ID(env->p[1]) = 1;
    MASS(env->p[1]) = 0.05;
    SOFT(env->p[1]) = DEFAULT_SOFTENING;

    env->p[0]->r[0] = -1; env->p[0]->r[1] = 0; env->p[0]->r[2] = 0;
    env->p[1]->r[0] = 1; env->p[1]->r[1] = 0; env->p[1]->r[2] = 0;

    double mass=0;
    cm[0] = cm[1] = cm[2] = 0;
    for (i=0; i < env->N; i++)
    {
        mass += MASS(env->p[i]);
        for (j=0; j < 3; j++) cm[j] += MASS(env->p[i]) * env->p[i]->r[j];
    }
    for (j=0; j < 3; j++) cm[j] /= mass;

    fprintf(err, "cm = (%f %f %f)\n", cm[0], cm[1], cm[2]);

    for (i=0; i < env->N; i++)
    {
        mass += MASS(env->p[i]);
        for (j=0; j < 3; j++) env->p[i]->r[j] -= cm[j];
    }

    env->p[0]->v[0] = 0; 
    env->p[0]->v[1] = sqrt(K2*MASS(env->p[1])*MAG(RADIUS(env->p[0])) / SQDIST(env->p[0], env->p[1]));
    env->p[0]->v[2] = 0;

    env->p[1]->v[0] = 0; 
    env->p[1]->v[1] = -env->p[0]->v[1] * (MASS(env->p[0]) / MASS(env->p[1]));
    env->p[1]->v[2] = 0;
}

void ic_threebody(struct env *env)
{
    int i;
    VL(1) fprintf(err, "Using three-body system.\n");

    env->N = 3;
    env->ps = (struct particle *)calloc(env->N, sizeof(struct particle));
    env->p  = (struct particle **)malloc(env->N * sizeof(struct particle *));

    for (i=0; i < 3; i++)
    {
        env->p[i] = &(env->ps[i]);
        PQ_INIT(env->p[i]);
        ID(env->p[i]) = i;
        MASS(env->p[i]) = 1;
        SOFT(env->p[i]) = DEFAULT_SOFTENING;
    }

    env->p[0]->r[0] = -0.995492; env->p[0]->r[1] = 0.0; env->p[0]->r[2] = 0.0;
    env->p[1]->r[0] =  0.995492; env->p[1]->r[1] = 0.0; env->p[1]->r[2] = 0.0;
    env->p[2]->r[0] =       0.0; env->p[2]->r[1] = 0.0; env->p[2]->r[2] = 0.0;

    env->p[0]->v[0] = -0.347902;
    env->p[0]->v[1] = -0.53393;
    env->p[0]->v[2] = 0;

    env->p[1]->v[0] = -0.347902;
    env->p[1]->v[1] = -0.53393;
    env->p[1]->v[2] = 0;

    env->p[2]->v[0] = 0.695804;
    env->p[2]->v[1] = 1.067860;
    //env->p[2]->v[1] = sqrt(K2*MASS(env->p[1])*MAG(RADIUS(env->p[2])) / SQDIST(env->p[2], env->p[1]));
    env->p[2]->v[2] = 0;

    for (i=0; i < env->N; i++)
        fprintf(err, "%i %f %f %f %f %f %f %f\n", 
            0, MASS(env->p[i]), 
            env->p[i]->r[0], env->p[i]->r[1], env->p[i]->r[2],
            env->p[i]->v[0], env->p[i]->v[1], env->p[i]->v[2]);
}

void ic_earthsun(struct env *env)
{
    if (verbosity > 0) fprintf(err, "Using earth/sun system.\n");

    env->N = 2;
    env->p = (struct particle **)malloc(env->N * sizeof(struct particle *));

    int S = 0;
    int E = 1;

    /* Sun */
    env->p[S] = (struct particle *)calloc(1, sizeof(struct particle));
    PQ_INIT(env->p[S]);
    ID(env->p[S]) = 1;
    MASS(env->p[S]) = 1 - 1e-6;
    SOFT(env->p[S]) = DEFAULT_SOFTENING;
    env->p[S]->r[0] = 0; env->p[S]->r[1] = 0; env->p[S]->r[2] = 0;
    env->p[S]->v[0] = 0; env->p[S]->v[1] = 0; env->p[S]->v[2] = 0;

    /* Earth */
    env->p[E] = (struct particle *)calloc(1, sizeof(struct particle));
    PQ_INIT(env->p[E]);
    ID(env->p[E]) = 0;
    MASS(env->p[E]) = 1e-6;
    SOFT(env->p[E]) = DEFAULT_SOFTENING;
    env->p[E]->r[0] = 1; env->p[E]->r[1] = 0; env->p[E]->r[2] = 0;
    env->p[E]->v[0] = 0; 
    env->p[E]->v[1] = -sqrt(K2*MASS(env->p[S]) / ABS(DIST(env->p[E], env->p[S]))); 
    env->p[E]->v[2] = 0;
}

void ic_randomsphere(struct env *env)
{
    int i;

    if (verbosity > 0) fprintf(err, "Using random system.\n");

    env->N = 2000;
    env->p = (struct particle **)malloc(env->N * sizeof(struct particle *));

    double rad = 10;
    double maxv = .0000000001;

    srand(time(NULL));
    for (i=0; i < env->N; i++)
    {
        env->p[i] = (struct particle *)calloc(1, sizeof(struct particle));
        PQ_INIT(env->p[i]);
        ID(env->p[i]) = i;
        MASS(env->p[i]) = 1;
        SOFT(env->p[i]) = DEFAULT_SOFTENING;

        env->p[i]->r[0] = rad * ((double)rand() / RAND_MAX) - rad/2;
        env->p[i]->r[1] = rad * ((double)rand() / RAND_MAX) - rad/2;
        env->p[i]->r[2] = rad * ((double)rand() / RAND_MAX) - rad/2;

        env->p[i]->v[0] = maxv * ((double)rand() / RAND_MAX) - maxv/2;
        env->p[i]->v[1] = maxv * ((double)rand() / RAND_MAX) - maxv/2;
        env->p[i]->v[2] = maxv * ((double)rand() / RAND_MAX) - maxv/2;
    }
}


#if 0
void write_state_every(double timeNext, int maxTime)
{
    int i;
    static int seqno=0;
    static int time=0;

    //qsort(env->p, env->N, sizeof(struct particle *), compar_radius);

    if ((time % maxTime) == 0)
    {
        time = 0;

        if (outfp[OUT_ENERGY])        write_energy(outfp[OUT_ENERGY], timeNext);
        if (outfp[OUT_DENSITY])       write_density(outfp[OUT_DENSITY]);
        if (outfp[OUT_MASS])          write_mass(outfp[OUT_MASS], timeNext);
        if (outfp[OUT_PARTICLE_INFO]) write_ascii_fp(outfp[OUT_PARTICLE_INFO]);
        if (outfp[OUT_TIMESTEenv->p])     write_timesteenv->p(outfp[OUT_TIMESTEenv->p]);

#ifdef ACCEL_HIST
        if (outfp[OUT_ACCEL_HIST])
        {
            //fseek(outfp[OUT_ACCEL_HIST], 0, SEEK_SET);
            for (i=0; i < ACCEL_HIST_SIZE; i++)
                fprintf(outfp[OUT_ACCEL_HIST], "%i %i\n", i-ACCEL_HIST_SIZE/2, accel_hist[i]);
            fprintf(outfp[OUT_ACCEL_HIST], "\n\n");
            fflush(outfp[OUT_ACCEL_HIST]);
        }
#endif

        if (save_format == SAVE_TIenv->p) write_tienv->p(outfilebase, seqno++);
        if (save_format == SAVE_ASCII) write_ascii(outfilebase, seqno++);

        if (verbosity > 0)
        {
#if 0
            if (ID(p) != oldid) 
            {
                oldid = ID(p);
                fprintf(err, "h(%i): %22.15e\n", ID(p), p->timeNext); //, h1, h2);
            }
#endif
            fprintf(err, "\r%22.15e %i   ", timeNext, AccelCount); //, h1, h2);

#if 0
            for (i=0; i < env->N; i++)
                fprintf(err, "%i ", env->p[i]->rung);
            fprintf(err, "\n");
#endif
        }
    }
}

void write_state(double timeNext)
{
    write_state_every(timeNext, 100);
}
#endif

