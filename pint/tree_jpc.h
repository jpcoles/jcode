#ifndef TREE_JS
#define TREE_JS

#include "pint.h"


struct tree 
{
    struct tree *left;
    struct tree *right;
    int d;
    real split;
    int count;

    char pad0[4];       
    /* 32 bytes */

    struct particle *list[MAX_PIC];
    /* 72 bytes */

    real max[3];
    real cr[3];

    real r[3];
    real mass;
    /* 40 bytes */

    char pad1[56];
};

void tree_build_jpc(struct env *env);
void tree_free_jpc(struct env *env);

#endif

