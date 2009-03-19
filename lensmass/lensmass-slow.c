#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <readline/readline.h>
#include "timing.h"

FILE *in  = stdin;
FILE *out = stdout;
FILE *err = stderr;

CPUDEFS

int main(int argc, char **argv)
{
    char *line;
    float **M = NULL;
    float *xs = NULL;
    int dim=0, old_dim=-1;
    int i, j, num_xs=0, num_models=0;

    float start_time=0, end_time=0, total_time=0;

    rl_instream = in;
    rl_outstream = out;

    /*========================================================================
     * Read commands or numerical input from the command line and update
     * the matrix accordingly, until there is no more input.
     *======================================================================*/
    while (1)
    {
        line = readline(NULL);
        if (line == NULL) break;
        if (line[0] == '#') continue;

        if (!strcmp("BEGIN NEW RUN", line))
        {
            if (M != NULL)
            {
                for (i=0; i < dim; i++) free(M[i]);
                free(M);
                M = NULL;
            }

            if (xs != NULL) 
            {
                free(xs);
                xs = NULL;
            }

            old_dim = -1;
        }
        if (!strcmp("BEGIN MODEL", line))
        {
            line = readline(NULL);
            if (line == NULL) break;

            dim = atoi(line);

            if (dim != old_dim)
            {
                if (M != NULL)
                {
                    for (i=0; i < dim; i++) free(M[i]);
                    free(M);
                    M = NULL;
                }

                if (xs != NULL) 
                {
                    free(xs);
                    xs = NULL;
                }

                old_dim = dim;
            }

            /*================================================================
             * Create a new matrix.
             *==============================================================*/
            if (M == NULL)
            {
                M = (float **)malloc(dim * sizeof(float *));
                for (i=0; i < dim; i++)
                {
                    M[i] = (float *)calloc(dim, sizeof(float));
                }
            }

            /*================================================================
             * Create an array for the x's.
             *==============================================================*/
            if (xs == NULL)
            {
                xs = (float *)malloc(dim * dim * sizeof(float));
            }

            num_xs = 0;

            start_time = CPUTIME;
        }
        else if (!strcmp("END MODEL", line))
        {
            num_models++;;

            /*================================================================
             * Now that the x's have been read, update the matrix.
             *==============================================================*/
            if (num_xs < dim)
                fprintf(err, "ERROR: Too few x's\n");
            else if (num_xs > dim)
                fprintf(err, "ERROR: Too many x's\n");
            else
            {
                /*================================================================
                 * Update the matrix such that the current state is
                 *      M_ij += xs[i] * xs[j];
                 *
                 * All the pointer mess actually does make a 2x speed up over
                 * indexing the array directly.
                 *==============================================================*/
                float **M_i = M;
                float *xs_i = xs;
                for (i=0; i < dim; i++, M_i++, xs_i++)
                {
                    float *M_ij = *M_i;
                    float *xs_j = xs;
                    register float xs_i_val = *xs_i;
                    for (j=0; j < dim; j++, xs_j++)
                    {
                        *M_ij++ += xs_i_val * *xs_j;
                    }
                }

                end_time = CPUTIME;
                total_time += end_time - start_time;

                fprintf(out, "[%ix%i matrix; %i models; %fs; %fs]\n", 
                        dim, dim, num_models, end_time - start_time, total_time);
            }
        }
        else if (xs != NULL)
        {
            xs[num_xs++] = atof(line);
        }

    }

    fprintf(out, "[%ix%i matrix; %i models]\n", dim, dim, num_models);
    for (i=0; i < dim; i++)
    {
        for (j=0; j < dim; j++)
        {
            fprintf(out, "%f ", M[i][j]);
        }
        fprintf(out, "\n");
    }

    return 0;
}

