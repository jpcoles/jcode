/*
 * 
 * Group id's begin with 1 (0 is for field particles)
 * Particle id's begin with 0
 *
 */
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <getopt.h>

const int debug_level = 0;
#define DBG(_lvl_) if (debug_level >= (_lvl_))

#define RELABEL 0

/* Halos must have at least n-1 parent halos to count. */
#define MIN_LEVEL -1

void help()
{
    fprintf(stderr, "Usage: ahf2tipgrp [--relabel] <#Particles> <AHF_particles> <grp-output>\n");
    exit(2);
}

int main(int argc, char **argv)
{
    int relabel = 0;
    int belong  = 0;
    size_t nParticles=0;
    FILE *in, *out;

    if (argc < 3) help();

    while (1)
    {
        int c;
        int option_index = 0;
        static struct option long_options[] = {
           {"relabel", 1, 0, 0},
           {"help", 1, 0, 'h'},
           {"belong", 1, 0, 'b'},
           {0, 0, 0, 0}
        };

        c = getopt_long(argc, argv, "h", long_options, &option_index);
        if (c == -1)
           break;

        switch (c) {
            case 0:
                if (!strcmp("relabel", long_options[option_index].name))
                    relabel = 1;
            case 'h':
                help();
                break;
            case 'b':
                belong = 1;
                break;
        }
    }

    if (argc-optind < 3) help();

    if ((nParticles = atol(argv[optind])) < 0)
    {
        fprintf(stderr, "Bad number of particles %s\n", argv[optind]);
        exit(2);
    }
    optind++;

    if (!strcmp("-", argv[optind]))
    {
        in = stdin;
    }
    else if ((in = fopen(argv[optind], "r")) == NULL)
    {
        fprintf(stderr, "Can't open input file %s\n", argv[optind]);
        exit(2);
    }
    optind++;

    if (!strcmp("-", argv[optind]))
    {
        out = stdout;
    }
    if ((out = fopen(argv[optind], "w")) == NULL)
    {
        fprintf(stderr, "Can't open output file %s\n", argv[optind]);
        exit(2);
    }
    optind++;

    int *list       = NULL;
    int *levelCount = NULL;

    size_t groupCountLen  = 0;
    int *groupCountBefore = NULL;
    int *groupCountAfter  = NULL;

    list       = (int *)calloc(nParticles, sizeof(int));
    levelCount = (int *)calloc(nParticles, sizeof(int));

    fprintf(stderr, "Reading input...\n");

    int i, gid=0, pid, nGrpParticles;
    while (!feof(in))
    {
        if (fscanf(in, "%i", &nGrpParticles) == 1)
        {
            if (gid == groupCountLen)
            {
                if (groupCountLen == 0) groupCountLen = 512;
                else groupCountLen *= 2;

                groupCountBefore = realloc(groupCountBefore, (groupCountLen+1) * sizeof(int));
                memset(groupCountBefore + gid, 0, (groupCountLen+1-gid) * sizeof(int));

                groupCountAfter = realloc(groupCountAfter, (groupCountLen+1) * sizeof(int));
                memset(groupCountAfter + gid, 0, (groupCountLen+1-gid) * sizeof(int));
            }

            gid++;

            if (nGrpParticles == 0)
            {
                fprintf(stderr, "Corrupt file? Group %i has 0 particles.\n", gid);
                exit(1);
            }

            for (i=0; i < nGrpParticles; i++)
            {
                if (fscanf(in, "%i", &pid) == 1)
                {
                    if (! (0 <= pid && pid < nParticles) )
                    {
                        fprintf(stderr, "ERROR: Particle has id %i.\n", pid);
                        exit(1);
                    }

#if 0
                    if (pid >= listLen) /* '>' because ids are not monotonically increasing */
                    {
                        if (listLen == 0) listLen = 512;
                        while (listLen <= pid) listLen *= 2;

                        list       = realloc(list,       listLen * sizeof(int));
                        levelCount = realloc(levelCount, listLen * sizeof(int));
                        assert(list != NULL);
                        assert(levelCount != NULL);

                        memset(      list+nParticles, 0, listLen-nParticles);
                        memset(levelCount+nParticles, 0, listLen-nParticles);
                    }

                    if (pid > nParticles) nParticles = pid;
#endif

#if 1
                    if (list[pid] != 0) 
                        groupCountAfter[list[pid]-1]--;

                    groupCountBefore[gid-1]++;
                    groupCountAfter[gid-1]++;

#if 0
                    if (gid == 1024)
                        fprintf(stderr, "** %i %i\n", groupCountBefore[gid-1], groupCountAfter[gid-1]);
#endif

                    list[pid] = gid;
#else
                    if (gid >= 3) list[pid] = gid;
                    else          list[pid] = 0;
#endif

                    levelCount[pid]++;
                }
            }

        }
    }

    if (gid == 0)
    {
        fprintf(stderr, "WARNING: No groups found!\n");
        //exit(1);
    }

    for (i=0; i < nParticles; i++)
    {
        if (list[i] > nParticles)
        {
            fprintf(stderr, "ERROR: Particle %i has id %i, but there are only %ld particles.\n",
                    i+1, list[i], nParticles);
            exit(1);
        }
    }

    for (i=0; i < nParticles; i++)
        if (levelCount[i] < MIN_LEVEL) list[i] = 0;

    fprintf(stderr, "Counting group membership...\n");
    int *groups = (int *)calloc((size_t)(gid+1), sizeof(int));  assert(groups != NULL);
    for (i=0; i < nParticles; i++)
        groups[list[i]]++;


    fprintf(stderr, "%i Groups.\n", gid);

    if (relabel)
    {
        fprintf(stderr, "Relabeling...\n");
        /* 
         * Some funny groups can lose all their particles to subgroups (how?). In 
         * this case, we must remove those groups and relabel all the others that
         * follow it in sequence. 
         */
        int *removed = (int *)calloc((size_t)(gid+1), sizeof(int));  assert(removed != NULL);

        int emptyGroups=0;
        for (i=1; i <= gid; i++)
        {
            emptyGroups += groups[i] == 0;
            removed[i] = emptyGroups;
        }

        for (i=0; i < nParticles; i++)
        {
            assert(groups[list[i]] != 0); /* Particles should always be somewhere */
            list[i] -= removed[list[i]];
        }

        fprintf(stderr, "WARNING: %i groups were empty and have been removed.\n", emptyGroups);
        fprintf(stderr, "WARNING: Group file and stat file are no longer consistent.\n");

        gid -= emptyGroups;

        fprintf(stderr, "%i Groups.\n", gid);
    }
    else
    {

        DBG(1)
        {
            memset(groups, 0, sizeof(int) * (gid+1));
            for (i=0; i < nParticles; i++)
                groups[list[i]]++;

            int sum=0;
            int sum2=0;
            for (i=0; i <= gid; i++)
            {
                fprintf(stderr, "[%4i] -> %i (%i %i)\n", i, groups[i], groupCountBefore[i-1], groupCountAfter[i-1]);
                if (i>1) sum += groupCountAfter[i-1];
                sum2 += groups[i];
            }

            fprintf(stderr, "nParticles: % 9ld\n", nParticles);
            fprintf(stderr, "Field:      % 9i\n", groups[0]);
            fprintf(stderr, "Main halo:  % 9i\n", groups[1]);
            fprintf(stderr, "Remainder:  % 9i\n", sum);
        }

        for (i=1; i <= gid; i++)
        {
            if (groups[i] == 0)
                fprintf(stderr, "Group %i has no particles left.\n", i);
        }
    }

    free(groups);

    fprintf(stderr, "Writing output...\n");

    fprintf(out, "%ld\n", nParticles);
    for (i=0; i < nParticles; i++)
        fprintf(out, "%i\n", list[i]);

    if (in  != stdin)  fclose(in);
    if (out != stdout) fclose(out);

    return 0;
}

