#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include <assert.h>
#include "pint.h"
#include "io.h"
#include "tipsy.h"
#include "leimkuhler.h"
#include "lk-tree.h"
#include "timing.h"

static CPUDEFS

#define TIME_CAP 2


extern int AccelCount;

extern int comparisonCount;
extern int nodeCount;

static double SELECT_TIME_STEP(struct pqNode *p)
{
    double eta = DFLT_ETA;
    double v = MAG(VEL(p));
    double T = 0.5*MASS(p)*v*v;
    return 0.25 * eta*eta / (T - p->U);
}

/* Update particle position */
static inline void DRIFT(struct pqNode *p, double dt)
{
    int i; _3x_(i) p->r[i] += p->v[i] * dt;
}

/* Update particle velocity */
static inline void KICK(struct pqNode *p, double dt)
{
    int i; _3x_(i) p->v[i] += p->a[i] * dt;
}

#if DYNAMIC_TREE
int walk_tree(pTree *tree, struct pqNode **ps, int index)
{
    if (tree == NULL) return index;

    fprintf(err, "( ");
    index = walk_tree(tree->pLeft, ps, index);

    struct pqNode *l;
    for (l = tree->cell.list; l != NULL; l = l->neighbor)
    {
        //ps[index++] = l;
        fprintf(err, "%.2f ", l->r[0]);
    }

    index = walk_tree(tree->pRight, ps, index);
    fprintf(err, ") ");

    return index;
}
#else
int walk_tree(pTree *tree, struct pqNode **ps, int index)
{
    if (tree == NULL) return index;

    fprintf(err, "( ");
    index = walk_tree(tree->pLeft, ps, index);

    if (tree->pLeft == NULL && tree->pRight == NULL)
    {
        int i;
        for (i=tree->iLower; i <= tree->iUpper; i++)
            fprintf(err, "%.2f ", ps[i]->r[0]);
    }

    index = walk_tree(tree->pRight, ps, index);

    fprintf(err, ") \n");

    return index;
}
#endif

static void lkACCEL_AND_POTENTIAL(struct pqNode *p, double time)
{
    int i, j;
    FLOAT a[3] = {0.0, 0.0, 0.0};
    FLOAT dr[3];

    AccelCount++;

    double potential = 0;

    /* Sum the forces from other particles */
    for (i=0; i < NPARTICLES; i++)
    {
        if (ID(p) == ID(ps[i])) continue;

        for (j=0; j < 3; j++)
            dr[j] = ps[i]->r[j] - p->r[j];

        //double dist  = 1./sqrt(SQMAG(dr) + pow(SOFTENING(p,ps[i]), 2));
        double dist  = 1. / sqrt(SQMAG(dr));

        double f = K2 * MASS(ps[i]);

        potential -= f * dist;

        for (j=0; j < 3; j++) 
        {
            double add = (((f * dr[j]) * dist) * dist) * dist;
            a[j] += add;
        }
    }

    /* Set acceleration */
    for (i=0; i < 3; i++) p->a[i] = a[i];

    p->U = potential;
}

void leimkuhler_tree(double maxtime)
{
    int i;
    const double tinterval = 1.0;
    float time_start, time_end;

    if (verbosity > 0)
    {
        fprintf(err, "Using leimkuhler_tree integrator.\n"
                     "TIME_CAP:      %3i\n", 
                     TIME_CAP);
        fprintf(err, "-------------------------------------------------\n");
    }

    /*========================================================================
     * Calculate the initial forces and time steps, start the integrator
     * by first drifting.
     *======================================================================*/
    for (i=0; i < NPARTICLES; i++)
        lkACCEL_AND_POTENTIAL(ps[i], 0);

    for (i=0; i < NPARTICLES; i++)
    {
        double h2 = sqrt(SELECT_TIME_STEP(ps[i]));
        DRIFT(ps[i], h2);
        ps[i]->time = 0;
        ps[i]->timeNext = h2;
    }
    write_state(0);

    /*========================================================================
     * Put each particle into the priority queue
     *======================================================================*/
    struct pqNode *head = ps[0];
    for (i=1; i < NPARTICLES; i++)
    {
        struct pqNode *p = ps[i];
        PQ_MERGE(head, p);
    }

    pTree *tree;

    tree = pintBuildTree(ps, NULL);

#if 1
    time_start = CPUTIME;
    for (i=0; i < 50000; i++)
        tree = pintBuildTree(ps, NULL);
    time_end   = CPUTIME;
#endif
#if DYNAMIC_TREE
    //fprintf(err, "tree->cell.fMass = %e\n", tree->cell.fMass);
#endif
    fprintf(err, "comparisonCount = %i\n", comparisonCount);
    fprintf(err, "nodeCount = %i\n", nodeCount);
    fprintf(err, "tree construction time = %e\n", (time_end - time_start) / 50000.);
    //fprintf(err, "tree construction time = %e\n", (time_end - time_start) );
    //assert(tree->cell.fMass == 1.0);

#if 0
#if DYNAMIC_TREE
    int len = walk_tree(tree, NULL, 0);
    fprintf(err, "\n");
#else
    int len = walk_tree(tree, ps, 0);
    fprintf(err, "\n");
#endif
#endif

#if 0
#if DYNAMIC_TREE
    fprintf(err, "\n\n");
    struct pqNode **ps2 = (struct pqNode **)malloc(NPARTICLES * sizeof(struct pqNode *));
    int len = walk_tree(tree, ps2, 0);
    fprintf(err, "len = %i\n", len);
    assert(len == NPARTICLES);

    struct pqNode *ps3 = (struct pqNode *)malloc(NPARTICLES * sizeof(struct pqNode));
    assert(ps3 != NULL);
    for (i=0; i < NPARTICLES; i++)
        memcpy(&(ps3[i]), ps2[i], sizeof(struct pqNode));

    for (i=0; i < NPARTICLES; i++)
        ps2[i] = &(ps3[i]);

    time_start = CPUTIME;
    for (i=0; i < 5000; i++)
        tree = pintBuildTree(ps2, NULL);
    time_end   = CPUTIME;
    fprintf(err, "tree->cell.fMass = %e\n", tree->cell.fMass);
    fprintf(err, "tree construction time = %e\n", (time_end - time_start) / 5000.);
#endif
#endif
    //assert(tree->cell.fMass == 1.0);

    return;

    /*========================================================================
     * Run the simulation
     *======================================================================*/
    double prevtimeNext = 0;
    double tout = tinterval;

    fflush(out);
    fflush(err);

    RESET_ACCEL_HIST();
    while (1)
    {
        struct pqNode *p;
        PQ_REMOVE_MIN(head, p);

        assert(prevtimeNext <= p->timeNext);
        prevtimeNext = p->timeNext;

        fprintf(out, "%i %f %f %f\n", p->id, p->r[0], p->r[1], p->r[2]);
        fflush(out);

        /*====================================================================
         * Write output
         *==================================================================*/
        if (p->timeNext < tout && tout <= head->timeNext) 
        {
            tout += tinterval;
            RESET_ACCEL_HIST();
            for (i=0; i < NPARTICLES; i++) ACCEL_HIST(ps[i]->a);

            write_state(p->timeNext);
            AccelCount = 0;

            if (p->timeNext < maxtime && maxtime <= head->timeNext) break;

#if !DYNAMIC_TREE
            //tree = pintBuildTree2(ps, tree);
            tree = pintBuildTree(ps, NULL);
#endif
        }

        /*====================================================================
         * Drift, Kick, Drift
         *==================================================================*/
        double h1 = p->timeNext - p->time;
        DRIFT(p, h1);
        lkACCEL_AND_POTENTIAL(p, p->timeNext);
        KICK(p, h1);
        double h2 = SELECT_TIME_STEP(p) / h1;
        KICK(p, h2);
        DRIFT(p, h2);

        /*====================================================================
         * Update particle's time and put it back in the queue
         *==================================================================*/
        p->time = p->timeNext + h2;
        p->timeNext = p->time + h2;

#if DYNAMIC_TREE
        _DA_ time_start = CPUTIME;
        tree = tree_update(p, tree);
        _DA_ time_end = CPUTIME;
        _DA_ fprintf(err, "tree update time = %e\n", time_end - time_start);
#endif

#if 0
#if DYNAMIC_TREE
        _DA_ time_start = CPUTIME;
        tree = tree_insert(p, tree);
        _DA_ time_end = CPUTIME;
        _DA_ fprintf(err, "tree insert time = %e\n", time_end - time_start);
#endif
#endif
        PQ_MERGE(head, p);
    }
}

