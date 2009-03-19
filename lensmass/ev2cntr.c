/*
 * Read the eigenvalue output from lensmass.c and sum the logs of the first
 * n values from each set. If the value for omega changes, add an additional
 * blank line, The output will look like 
 *
 *                  omega0 lambda sum0
 *                  omega0 lambda sum1
 *                  omega0 lambda sum2
 *
 *                  omega1 lambda sum3
 *                  omega1 lambda sum4
 *                  omega1 lambda sum5
 *                         ...
 *
 * This output can be given directly to gnuplot to plot contour maps.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <readline/readline.h>

FILE *in, *out, *err;

void help()
{
    fprintf(err, "ev2cntr <sum length>\n");
    exit(2);
}

int main(int argc, char **argv)
{
    char *line;
    int line_count;
    int ev_count=0;
    int read_ev = 0;
    float log_sum=0;
    int product_len=0;
    float omega=0, lambda=0;
    int dim = 0;

    in = stdin;
    out = stdout;
    err = stderr;

    if (argc < 2)
    {
        help();
    }

    product_len = atoi(argv[1]);

    line_count = 0;
    while (1)
    {
        line = readline(NULL);
        if (line == NULL) break;

        line_count++;

        if (strstr(line, "#BEGIN EIGENVALUES") == line)
        {
            float o=0;

            sscanf(line, "#BEGIN EIGENVALUES dim %i omega %f lambda %f", 
                &dim, &o, &lambda);

            if (o != omega)
                fprintf(out, "\n");

            omega = o;

            if (dim == 0)
                fprintf(out, "%f %f 0\n", omega, lambda);

            ev_count = 0;
            read_ev = 1;
            log_sum = 0;
        }
        else if (strstr(line, "#END EIGENVALUES") == line)
        {
            read_ev = 0;
        }
        else if (line[0] == '#') 
        {
            continue;
        }
        else if (read_ev)
        {
            if (ev_count < product_len)
            {
                ev_count++;
                double p = atof(line);

                log_sum += log(p);
                if (ev_count == product_len)
                {
                    fprintf(out, "%f %f %.60f\n", omega, lambda, log_sum);
                }
            }
        }
    }

    return 0;
}

