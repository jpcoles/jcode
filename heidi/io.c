#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "heidi.h"

Environment env;

int load_particles(const char *filename)
{
    int i, nTypes=1;
    FILE *fp;

    VL(2) fprintf(env.io.out, "Loading Particles from %s\n", filename);

    fp = fopen(filename, "r");

    char buf[256];

    char keyword[256], val[256];

    while (1)
    {
        fgets(buf, 256, fp);
        if (feof(fp)) break;
        if (buf[0] == 0 || buf[0] == '#') continue;

        if (sscanf(buf, "%s %s\n", keyword, val) < 2) continue;

        if (!strcmp("HeidiVersion:", keyword))
        {
        }
        else if (!strcmp("NumBodies:", keyword))
        {
            NP = strtol(val, NULL, 10);
        }
        else if (!strcmp("NumTypes:", keyword))
        {
            nTypes = strtol(val, NULL, 10);
        }
        else
            break;
    }

    env.plist.p = PTR_MALLOC(Particle, NP);
    env.particle_type_data = PTR_MALLOC(ParticleTypeData, nTypes + 1);

    env.particle_type_data[0].m = 0;
    env.particle_type_data[0].q = 0;

    for (i=0; i < NP; i++)
    {
        int   type;
        int   Z, N, E;
        float xx, xy, xz;
        float vx, vy, vz;
        float ax, ay, az;

        int nitems = 
        sscanf(buf, "%i %i %i %i %e %e %e %e %e %e %e %e %e",
            &type, 
            &Z, &N, &E,
            &xx, &xy, &xz, 
            &vx, &vy, &vz, 
            &ax, &ay, &az);

        if (nitems != 13) break;

        P_type(i) = type;

        if (type != 0 && !(Z == 0 && N == 0 && E == 0))
        {
            P_m(i) = (Z * _Mp) + (N * _Mn) + (E * _Me);
            P_q(i) = (Z * _Qp) + (N * _Qn) + (E * _Qe);
            P_Z(i) = Z;
            P_N(i) = N;
            P_E(i) = E;
        }

        P_xx(i) = xx;
        P_xy(i) = xy;
        P_xz(i) = xz;

        P_vx(i) = vx;
        P_vy(i) = vy;
        P_vz(i) = vz;

        P_ax(i) = ax;
        P_ay(i) = ay;
        P_az(i) = az;

        fgets(buf, 256, fp);
        if (feof(fp)) break;
        if (buf[0] == 0 || buf[0] == '#') continue;
    }

    fclose(fp);

    if (i != NP)
    {
        VL(1) fprintf(env.io.err, 
                      "ERROR: Expected %i particles, but only read %i.\n",
                      NP, i);

        return 1;
    }

    VL(2) 
    {
        fprintf(stderr, "NumBodies: %i\n", NP);
        fprintf(stderr, "NumTypes: %i\n", nTypes);
        for (i=0; i <= nTypes; i++)
        {
            fprintf(stderr, "Type: %i Mass: %e Charge: %e\n", i, env.particle_type_data[i].m, env.particle_type_data[i].q);
        }
    }

    return 0;
}

void write_state()
{
    int i;

    char filename[sizeof(env.io.state_file) * 2];

    if (strlen(env.io.state_file) == 0) return;

    //static int pad = 5;
    
    //if (pad == 0) pad = (int)log10(env.t_end / env.t_step) + 1;

    if (env.io.state_fp == NULL)
    {
        snprintf(filename, sizeof(filename), "%s-state", env.io.state_file);

        env.io.state_fp = fopen(filename, "w");
        if (env.io.state_fp == NULL) return;
    }

    for (i=0; i < 1; i++)
    {
        //fprintf(fp, "%e %i %e %e %e %e %e %e %e %e %e\n", 
        fprintf(env.io.state_fp, 
            "%e %i %f %f %f %f %f %f %24.15e %24.15e %24.15e\n", 
            env.t_now, i, 
            P_xx(i), P_xy(i), P_xz(i),
            P_vx(i), P_vy(i), P_vz(i),
            P_ax(i), P_ay(i), P_az(i)
            );
    }

    fflush(env.io.stats_fp);

    if (env.io.stats_fp == NULL)
    {
        snprintf(filename, sizeof(filename), "%s-stats", env.io.state_file);

        env.io.stats_fp = fopen(filename, "w");
        if (env.io.stats_fp == NULL) return;
    }

    fprintf(env.io.stats_fp, 
        "%e %e \n",
        env.t_now, env.E.kin
        );

    fflush(env.io.state_fp);
}

int is_stdio(FILE *fp)
{
    return (fp == stdin) 
        || (fp == stdout)
        || (fp == stderr);
}

int io_close(FILE *fp)
{
    if (fp == NULL) return 0;
    if (is_stdio(fp)) return 0;

    return fclose(fp);
}

