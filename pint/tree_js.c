#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "pint.h"
#include "timing.h"
#include "tree_js.h"


#if USE_TREE_JS


static CPUDEFS

#define USE_MEMORY_POOL 1

#define NODE(n) (&(node_pool_ptr[n]))

/* the order here is important, as it is assumed elsewhere */
enum { SPLIT_X_DIM = 0, SPLIT_Y_DIM, SPLIT_Z_DIM }; 

/* Things for the memory allocation scheme */
#define NODE_POOL_SEGMENT_SIZE (10000)
static struct tree **node_pool = NULL;
static int nodePool_seg = 0; 
static int nodePool_off = NODE_POOL_SEGMENT_SIZE;

static double cellSize = 0;

static int comparisonCount = 0;
static int nodeCount=0;

static double bottom_up_time;
static int bottom_up_ncalls;

static double split_time;
static int split_ncalls;

static double top_down_time;
static int top_down_ncalls;

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

#define EXTEND_BOUNDS(R, MIN, MAX) \
{      if (R < MIN) MIN = R; \
  else if (R > MAX) MAX = R; }
static void crop(struct particle **ps, struct tree *node, int front, int back)
{
    int j;
    real min[3],max[3];

    //fprintf(err, "crop(): %i %i\n", front, back);

    for (j=0; j < 3; ++j) 
    {
        min[j] = ps[front]->r[j];
        max[j] = ps[front]->r[j];
    }

    while (++front <= back)
    {
        for (j=0; j < 3; ++j) 
            EXTEND_BOUNDS(ps[front]->r[j], min[j], max[j]);
    }

    for (j=0; j < 3; ++j) {
        node->bnd.r[j] = 0.5 * (max[j] + min[j]);
        node->bnd.max[j]    = 0.5 * (max[j] - min[j]);
    }

    //fprintf(err, "crop(): min=(%f %f %f) max=(%f %f %f) r=(%f %f %f)\n",
    //    min[0], min[1], min[2],
    //    max[0], max[1], max[2],
    //    node->bnd.r[0],
    //    node->bnd.r[1],
    //    node->bnd.r[2]);

}
#undef EXTEND_BOUNDS

static void split(struct tree *node, int *d, real *split)
{
    if (node->bnd.max[0] < node->bnd.max[1])
        *d = node->bnd.max[1] < node->bnd.max[2] ? SPLIT_Z_DIM : SPLIT_Y_DIM;
    else
        *d = node->bnd.max[0] < node->bnd.max[2] ? SPLIT_Z_DIM : SPLIT_X_DIM;

    *split = node->bnd.r[*d];

    //fprintf(err, "split(): %i %f\n", *d, *split);
}


#define NEW_NODE(parent, child, l, u) \
{\
    parent->child = alloc_node(); \
    parent->child->left = parent->child->right = NULL; \
    parent->child->iLower = l; \
    parent->child->iUpper = u; \
}

struct tree **sortStack = NULL;     /* this is the stack */

int partitionCount;

static int sort(struct env *env, int M)
{
    int sp=0;
    struct tree *left=NULL, *right=NULL;
    int nBucket = 0;

    struct particle **ps = env->p;
    struct tree *node = env->tree;

    /*========================================================================
     * Allocate stack
     *======================================================================*/
    int ns = (int)MMAX(1, floor(log(((double)(env->N+1))/(M+1))/log(2.0)));
    sp = 0;
    sortStack = (struct tree **)realloc(sortStack, ns * sizeof(struct tree *));
    assert(sortStack != NULL);

    if (node->iUpper - node->iLower + 1 <= M) 
        return 0;

    while (1)
    {
        int i, p, k;
        struct particle *t;
        int nl, nr;
        int d;
        real fSplit;

        assert(node != NULL);
        split(node, &d, &fSplit);

#if 0
        i = node->iLower;
        int j = node->iUpper;

        //partitionCount++;
        //PARTITION(ps, t, ->r[d], i,j, < fSplit,> fSplit);
    	while (i <= j && ((ps[i] ->r[d]) < fSplit)) { ++i; } 

#if 0
        while (i <= j)
        {
            switch ((i-j) & 0x7)
                case 7: if 

        }
#endif


    	while (i <= j && ((ps[j] ->r[d]) > fSplit)) { --j; } 
    	while (i < j) { 
            SWAP(ps[i], ps[j], t); 
            while ((ps[++i] ->r[d]) < fSplit) { } 
            while ((ps[--j] ->r[d]) > fSplit) { } 
    	}

#else

#       define PART(A, i, j, t, d, CMPL) \
        do { \
            int k = i--; \
            do { \
                if (A[k]->r[d] CMPL) { i++; SWAP(A[i], A[k], t); } \
            } while (++k <= j); \
            i++; \
        } while(0)

        i = node->iLower;
        const int j = node->iUpper;

        switch (d)
        {
            case 0: PART(ps, i, j, t, 0, < fSplit); break;
            case 1: PART(ps, i, j, t, 1, < fSplit); break;
            case 2: PART(ps, i, j, t, 2, < fSplit); break;
        }
#endif

        //fprintf(err, "(%i %i)\n", i, node->iUpper);

        nl = i - node->iLower;
        nr = node->iUpper - i + 1;

        //fprintf(err, "nl=%i nr=%i\n", nl, nr);

        /*========================================================================
         * If both sides of the partition are not empty then create two new nodes
         * and crop the bounding boxes. Otherwise, undo the split operation.
         *======================================================================*/
        if (nl > 0 || nr > 0) {  // jpc changed this from && to || so that a split is always created.
            NEW_NODE(node, left, node->iLower, i - 1);
            NEW_NODE(node, right, i, node->iUpper);
            nodeCount += 2;

            X3( node->left->bnd.r   [_] = node->bnd.r   [_] );
            X3( node->left->bnd.max [_] = node->bnd.max [_] );

            node->left->bnd.max[d]          *= 0.5;
            node->left->bnd.r[d]            -= node->left->bnd.max[d];
            left = node->left;

            X3( node->right->bnd.r   [_] = node->bnd.r   [_] );
            X3( node->right->bnd.max [_] = node->bnd.max [_] );

            node->right->bnd.max[d]         *= 0.5;
            node->right->bnd.r[d]           += node->right->bnd.max[d];

            right = node->right;
        }
        else
        {
            node->bnd.max[d] *= 0.5;
            if (nl > 0) {
                node->bnd.r[d] -= node->bnd.max[d];
                left = node;
            }
            else {
                node->bnd.r[d] += node->bnd.max[d];
                right = node;
            }
        }

        /*========================================================================
         * Now figure out which subfile to process next
         *======================================================================*/
        if (nl > M && nr > M) 
        {
            if (nr > nl) sortStack[sp++] = right, node = left;
            else         sortStack[sp++] = left,  node = right;
        }
        else 
        {
            if (nl > M) node = left;
            //else if (nl > 0) left->iLower = 0;

            if (nr > M) node = right;
            //else if (nr > 0) right->iLower = 0;
        }

        if (nl <= M && nr <= M) {
            if (sp) node = sortStack[--sp];      /* pop tn */
            else break;
        }

    }

    return 0;
}

static void populate(struct particle **ps, struct tree * const node)
{
    real m, mass, x, y, z;

    if (node == NULL) return;

    populate(ps, node->left);
    populate(ps, node->right);

    m = mass = x = y = z = 0;

    if (node->left == NULL  &&  node->right == NULL)
    {
        /*====================================================================
         * Calculate total mass and CoM for a leaf node
         *==================================================================*/
        int pj = node->iLower;
        //ps[pj].iBucket = iNode;
        for (pj = node->iLower; pj <= node->iUpper; ++pj) 
        {
            //p[pj].iBucket = iNode;
            m = ps[pj]->mass;
            mass += m;
            x += m * ps[pj]->r[0];
            y += m * ps[pj]->r[1];
            z += m * ps[pj]->r[2];
        }
        m = 1/mass;
        node->r[0]  = m * x;
        node->r[1]  = m * y;
        node->r[2]  = m * z;
        node->mass = mass;
    }
    else
    {
        /*====================================================================
         * Calculate total mass and CoM for a node
         *==================================================================*/
        struct tree *left  = node->left;
        struct tree *right = node->right;

        if (left != NULL)
        {
            mass = left->mass;
            x = left->mass * left->r[0];
            y = left->mass * left->r[1];
            z = left->mass * left->r[2];
        }

        if (right != NULL)
        {
            mass += right->mass;
            x += right->mass * right->r[0];
            y += right->mass * right->r[1];
            z += right->mass * right->r[2];
        }

        m = 1/mass;
        node->r[0]  = m * x;
        node->r[1]  = m * y;
        node->r[2]  = m * z;
        node->mass = mass;
    }
}

void tree_build_js(struct env *env) //struct particle **ps, struct tree *head)
{
    env->tree = NULL;
    if (env->tree == NULL)
    {
        env->tree = alloc_node();
        assert(env->tree != NULL);
        env->tree->iLower = 0;
        env->tree->iUpper = env->N - 1;
        env->tree->left = env->tree->right = NULL;
    }

    comparisonCount = 0;
    nodeCount = 1;

    crop(env->p, env->tree, 0, env->N - 1);

    sort(env, MAX_PIC);
    //populate(ps, head);

    VL(2) fprintf(err, "nodeCount = %i\n", nodeCount);
    VL(2) fprintf(err, "head->mass = %f\n", env->tree->mass);
    VL(2) fprintf(err, "head->r = %f %f %f\n", env->tree->r[0], env->tree->r[1], env->tree->r[2]);
}

void tree_free_js(struct env *env)
{
    int i;

    VL(1) fprintf(err, "tree_free_js: Deallocating tree.\n");
    if (sortStack != NULL) { free(sortStack); sortStack = NULL; }
    if (node_pool != NULL)
    {
        for (i=0; i < nodePool_seg; i++)
            if (node_pool[i] != NULL) free(node_pool[i]);
        free(node_pool);
        node_pool = NULL;
    }
    
    nodePool_seg = 0; 
    nodePool_off = NODE_POOL_SEGMENT_SIZE;

    env->tree = NULL;
}

#endif
