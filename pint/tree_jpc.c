#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "pint.h"
#include "timing.h"
#include "tree_jpc.h"

#if USE_TREE_JPC

static CPUDEFS

#define USE_MEMORY_POOL 1

enum { SPLIT_X_DIM = 0, SPLIT_Y_DIM, SPLIT_Z_DIM }; /* the order here is important, as it is assumed elsewhere */

/* Things for the memory allocation scheme */
#define NODE_POOL_SEGMENT_SIZE (10000)
static struct tree **node_pool = NULL;
static int nodePool_seg = 0; 
static int nodePool_off = NODE_POOL_SEGMENT_SIZE;

double cellSize = 0;

int comparisonCount = 0;
int nodeCount=0;

double bottom_up_time;
int bottom_up_ncalls;

double split_time;
int split_ncalls;

double top_down_time;
int top_down_ncalls;

#define FABS(x) (((x) < 0) ? -(x) : (x))

/*============================================================================
 * Return a pointer to an empty (zero'd) struct tree.
 *
 * The node pool is a collection of pointers to arrays of struct tree's. There 
 * are three reasons why we don't simply realloc() a single large array. First,
 * doing so may not be feasible in the case where a contiguous length can
 * not be found. If space can be found, moving a data from one region to
 * another is expensive and not necessary in our case. Third, and perhaps
 * the most important, the pointers returned from this function are stored in
 * other parts of the code. If the array changes location then the values of
 * the pointers change and the rest of the code has stale pointers.
 *==========================================================================*/
static struct tree *alloc_node()
{
    struct tree *ptr;
    static struct tree *node_pool_ptr;

    if (nodePool_off == NODE_POOL_SEGMENT_SIZE)
    {
        nodePool_off = 0;

        node_pool
            = (struct tree **)realloc(node_pool, (nodePool_seg+1) * sizeof(struct tree *));

        node_pool[nodePool_seg]
            = node_pool_ptr
            = (struct tree *)malloc(NODE_POOL_SEGMENT_SIZE * sizeof(struct tree));

        assert(node_pool != NULL);

        nodePool_seg++;
    }

    nodePool_off++;

    ptr = node_pool_ptr++;
    assert(ptr != NULL);
    memset(ptr, 0, sizeof(struct tree));

    return ptr;
}


static void free_node(struct tree *ptr)
{
#if !USE_MEMORY_POOL
    free(ptr);
#endif
}

#if 0
static void split2(struct tree *node, int *d, real *split)
{
    if (node->max[0] < node->max[1])
        *d = node->max[1] < node->max[2] ? SPLIT_Z_DIM : SPLIT_Y_DIM;
    else
        *d = node->max[0] < node->max[2] ? SPLIT_Z_DIM : SPLIT_X_DIM;

    *split = node->cr[*d];

    //fprintf(err, "split(): %i %f\n", *d, *split);
}

double inbounds_time;
int inbounds_ncalls;

int inbounds(struct particle *p, struct tree *node)
{
    //_DA_ fprintf(err, "inbounds(): %f %f %f\n", p->r[0], p->r[1], p->r[2]);
    //_DA_ fprintf(err, "            %f %f %f\n", cell->r[0], cell->r[1], cell->r[2]);
    //_DA_ fprintf(err, "            %f %f %f\n", cell->max[0], cell->max[1], cell->max[2]);

    //inbounds_ncalls++;

    return (FABS(p->r[0] - node->cr[0]) <= node->max[0] + 1e-5) 
        && (FABS(p->r[1] - node->cr[1]) <= node->max[1] + 1e-5) 
        && (FABS(p->r[2] - node->cr[2]) <= node->max[2] + 1e-5);
}
#endif

static void assign_split(struct tree *node, int *d0, real *split0)
{
    int d;

    if (node->max[0] < node->max[1])
        d = node->max[1] < node->max[2] ? SPLIT_Z_DIM : SPLIT_Y_DIM;
    else
        d = node->max[0] < node->max[2] ? SPLIT_Z_DIM : SPLIT_X_DIM;

    //assert(d != 3);
    node->split = node->cr[d];
    node->d = d;

    *d0 = d;
    *split0 = node->cr[d];
}

#if 0
struct particle *join_lists(struct particle *a, struct particle *b)
{
    struct particle *prev, *t;

    if (a == NULL) return b;
    if (b == NULL) return a;

    for (prev = t = a; t != NULL; prev = t, t = t->neighbor) {}

    prev->neighbor = b;

    return a;
}
#endif

#define CREATE_BRANCH(head, d, T, OP) \
{\
    head-> T  = alloc_node();    \
    head-> T -> d = 4; \
    X3( head-> T ->cr[_]   = head->cr[_] );  \
    X3( head-> T ->max[_]  = head->max[_] ); \
    head-> T ->max[d]     *= 0.5F;                       \
    head-> T ->cr[d]       = head-> T ->cr[d] OP head-> T ->max[d];   \
    nodeCount++; \
}
    //assign_split(head-> T );                                   


#if 0
struct tree *split_cell(const struct particle *p, struct tree *node)
{
    int i, d;
    real split;

    for (;;)
    {
        assign_split(node, &d, &split);

        CREATE_BRANCH(node, d, left, -);
        CREATE_BRANCH(node, d, right, +);

        struct tree *top = node;
        struct tree *left = node->left;
        struct tree *right = node->right;

        int nl=0, nr=0;
        for (i=0; top->list[i] != NULL; i++)
        {
            if (top->list[i]->r[d] < split)
                left->list[nl++] = top->list[i];
            else
                right->list[nr++] = top->list[i];
        }

        left->list[nl] = NULL;
        right->list[nr] = NULL;
        left->count  = nl;
        right->count = nr;

        node->list[0] = NULL;
        node->count = 0;

        if (p->r[d] < split)
        {
            node = node->left;
            if (nr != 0) break;
        }
        else
        {
            node = node->right;
            if (nl != 0) break;
        }

    }

    return node;
}
#endif

void tree_build(struct env *env)
{
    VL(1) fprintf(err, "Finding leaf...\n");
    int i;

    struct particle *p = env->ps;
    for (i=env->N - 1; i >= 0; i--, p++)
    {
        struct tree *node = env->tree;

        while (node->left != NULL) // && node->right != NULL)
        {
            register float t = p->r[node->d] - node->split;
            int q0 = 0-signbit(t);
            int q1 = signbit(t)-1;

            long r0 = (long)(node->left);
            long r1 = (long)(node->right);
            long q2 = r0 & q0;
            long q3 = r1 & q1;
            node = (struct tree *)(q2 + q3);

            //node = (struct tree *)(((unsigned long)(node->left)  &  q1)
             //                    + ((unsigned long)(node->right) & ~q1));
        }

        while (node->count == MAX_PIC) 
        {
            int d, j;
            real split;

            assign_split(node, &d, &split);

            CREATE_BRANCH(node, d, left, -);
            CREATE_BRANCH(node, d, right, +);

            struct tree *top = node;
            struct tree *left = node->left;
            struct tree *right = node->right;
            struct particle **tlist = top->list;

            int nl=0, nr=0;
            //for ( ; tlist != NULL; tlist++)

#define SELECT_D(d) \
switch (d) { \
    case 0: SELECT_COUNT(0); break; \
    case 1: SELECT_COUNT(1); break; \
    case 2: SELECT_COUNT(2); break; }

#define SELECT_COUNT(y) \
switch (top->count) { \
    case 8: GO_LEFT_OR_RIGHT(7, y); /* no break */ \
    case 7: GO_LEFT_OR_RIGHT(6, y); /* no break */ \
    case 6: GO_LEFT_OR_RIGHT(5, y); /* no break */ \
    case 5: GO_LEFT_OR_RIGHT(4, y); /* no break */ \
    case 4: GO_LEFT_OR_RIGHT(3, y); /* no break */ \
    case 3: GO_LEFT_OR_RIGHT(2, y); /* no break */ \
    case 2: GO_LEFT_OR_RIGHT(1, y); /* no break */ \
    case 1: GO_LEFT_OR_RIGHT(0, y); /* no break */ }

#define GO_LEFT_OR_RIGHT(x,y) { \
if (tlist[ x ]->r[ y ] < split) \
    left->list[nl++] = tlist[ x ]; \
else \
    right->list[nr++] = tlist[ x ]; }
                
            j = top->count;
            SELECT_D(d);

            left->count  = nl;
            right->count = nr;

            top->count = 0;
            top->list[0] = NULL;

            register unsigned long q1 = -(p->r[d] < split);

            node = (struct tree *)(((unsigned long)(node->left)  &  q1)
                                 + ((unsigned long)(node->right) & ~q1));

        }

        node->list[node->count++] = p;
    }

}

#define EXTEND_BOUNDS(R, MIN, MAX) \
{      if (R < MIN) MIN = R; \
  else if (R > MAX) MAX = R; }
static void crop(struct particle **ps, struct tree * node, int front, int back)
{
    real min[3],max[3];

    //fprintf(err, "crop(): %i %i\n", front, back);

    X3( min[_] = ps[front]->r[_] );
    X3( max[_] = ps[front]->r[_] );

    while (++front <= back)
        X3( EXTEND_BOUNDS(ps[front]->r[_], min[_], max[_]) );

    X3( node->cr[_]  = 0.5 * (max[_] + min[_]) );
    X3( node->max[_] = 0.5 * (max[_] - min[_]) );
        //assert(node->max[j] >= 0);

#if 0
    fprintf(err, "crop(): min=(%f %f %f) max=(%f %f %f) fCenter=(%f %f %f)\n",
        min[0], min[1], min[2],
        max[0], max[1], max[2],
        node->cr[0],
        node->cr[1],
        node->cr[2]);
#endif

}
#undef EXTEND_BOUNDS

void tree_build_jpc(struct env *env)
{
    int i=0, d;
    real split;

    if (env->tree == NULL)
    {
        env->tree = alloc_node();
        env->tree->d = 4;
    }

    crop(env->p, env->tree, 0, env->N-1);
    assign_split(env->tree, &d, &split);

    comparisonCount = 0;
    nodeCount = 1;

    fprintf(err, "sizeof(struct particle) = %ld\n", sizeof(struct particle));
    fprintf(err, "sizeof(struct tree) = %ld\n", sizeof(struct tree));
    fprintf(err, "sizeof(struct particle *) = %ld\n", sizeof(struct particle *));
    tree_build(env);

    fprintf(err, "nodeCount = %i\n", nodeCount);
}

void tree_free_jpc(struct env *env)
{
    VL(1) fprintf(err, "tree_free_jps: Deallocating tree.\n");
    if (node_pool != NULL)
    {
        int i;
        for (i=0; i < nodePool_seg; i++)
            if (node_pool[i] != NULL) free(node_pool[i]);
        free(node_pool);
        node_pool = NULL;
    }
    
    nodePool_seg = 0; 
    nodePool_off = NODE_POOL_SEGMENT_SIZE;

    env->tree = NULL;
}

#if 0
    int i;
    struct tree t;

    double start = CPUTIME, end;

    t.bnd.max[0] = 1;
    t.bnd.max[1] = 2;
    t.bnd.max[2] = 3;
    for (i=0; i < 10000000; i++) assign_split(&t);

    t.bnd.max[0] = 2;
    t.bnd.max[1] = 1;
    t.bnd.max[2] = 3;
    for (i=0; i < 10000000; i++) assign_split(&t);

    t.bnd.max[0] = 3;
    t.bnd.max[1] = 2;
    t.bnd.max[2] = 1;
    for (i=0; i < 10000000; i++) assign_split(&t);

    t.bnd.max[0] = 3;
    t.bnd.max[1] = 1;
    t.bnd.max[2] = 2;
    for (i=0; i < 10000000; i++) assign_split(&t);

    end = CPUTIME;

    fprintf(out, "time: %f\n", (end - start));

}
#endif


#endif
