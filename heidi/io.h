#ifndef IO_H
#define IO_H

#include <stdio.h>

typedef struct
{
    FILE *in, *out, *err;
    FILE *state_fp;
    FILE *stats_fp;

    int verbosity;

    char state_file[256];
} IO;

int load_particles(const char *filename);
void write_state();
int io_close(FILE *fp);
int is_stdio(FILE *fp);

#endif
