#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "pint.h"
#include "lk-tree.h"
#include "timing.h"

static CPUDEFS

#define USE_MEMORY_POOL 1

#define NODE(n) (&(node_pool_ptr[n]))

enum { SPLIT_X_DIM = 0, SPLIT_Y_DIM, SPLIT_Z_DIM }; /* the order here is important, as it is assumed elsewhere */

pTree * node_pool = NULL;
pTree * node_pool_ptr = NULL;

long nodePoolSize=0;
long nodePoolUsed=0;

struct pList_s *item_pool = NULL;
int itemPoolSize=0;
int itemPoolUsed=0;

double cellSize = 0;

int comparisonCount = 0;
int nodeCount=0;

double bottom_up_time;
int bottom_up_ncalls;

double split_time;
int split_ncalls;

double top_down_time;
int top_down_ncalls;


#if 0
double inline FABS(double a)
{
    if (a < 0) return -a;
    return a;
}
#else
//#define FABS(x) ((FLOAT)(*(&(*((int)(&x)) & ~(1 << (sizeof(FLOAT)-1))))))
#define FABS(x) (((x) < 0) ? -(x) : (x))
#endif

/* Return a pointer to an empty (zero'd) pTree */
static TREE_TYPE alloc_node()
{
    TREE_TYPE ptr;


#if USE_MEMORY_POOL
    if (nodePoolUsed == nodePoolSize)
    {
        fprintf(err, "\n--------------------------------\n");
        fprintf(err, "sizeof(void *) = %ld\n", sizeof(void *));
        fprintf(err, "sizeof(FLOAT) = %ld\n", sizeof(FLOAT));
        fprintf(err, "sizeof(unsigned int) = %ld\n", sizeof(unsigned int));
#if DYNAMIC_TREE
        fprintf(err, "sizeof(union split_u) = %ld\n", sizeof(union split_u));
        fprintf(err, "sizeof(cell_t) = %ld\n", sizeof(cell_t));
#endif
        fprintf(err, "sizeof(pTree) = %ld\n", sizeof(pTree));
        nodePoolSize += (int)log(NPARTICLES) * 150;
        node_pool = (pTree *)realloc(node_pool, nodePoolSize * sizeof(pTree) + 64);
        node_pool_ptr = (pTree *)((unsigned long)node_pool + 64 - ((unsigned long)node_pool & 63));

        fprintf(err, "Allocated %i pTree objects (%i + 64 bytes)\n", nodePoolSize, nodePoolSize * sizeof(pTree));
        fprintf(err, "node_pool @ %p\n", node_pool);
        fprintf(err, "node_pool_ptr @ %p\n", node_pool_ptr);

        memset(node_pool_ptr + nodePoolUsed, 0, (nodePoolSize - nodePoolUsed) * sizeof(pTree));
        assert(node_pool != NULL);
        fprintf(err, "--------------------------------\n");
    }

    //ptr = node_pool_ptr + nodePoolUsed;
    ptr = nodePoolUsed;
    nodePoolUsed++;

    //memset(ptr, 0, sizeof(pTree));
#else
    ptr = (pTree *)calloc(1, sizeof(pTree));
    assert(ptr != NULL);
#endif

    return ptr;
}

static void free_node(pTree *ptr)
{
#if !USE_MEMORY_POOL
    free(ptr);
#endif
}

#if !DYNAMIC_TREE

#define EXTEND_BOUNDS(R, MIN, MAX) \
{      if (R < MIN) MIN = R; \
  else if (R > MAX) MAX = R; }
static void crop(struct pqNode **ps, pTree *node, int front, int back)
{
    int j;
    FLOAT fMin[3],fMax[3];

    //fprintf(err, "crop(): %i %i\n", front, back);

    for (j=0; j < 3; ++j) 
    {
        fMin[j] = ps[front]->r[j];
        fMax[j] = ps[front]->r[j];
    }

    while (++front <= back)
    {
        for (j=0; j < 3; ++j) 
            EXTEND_BOUNDS(ps[front]->r[j], fMin[j], fMax[j]);
    }

    for (j=0; j < 3; ++j) {
        node->bnd.fCenter[j] = 0.5 * (fMax[j] + fMin[j]);
        node->bnd.fMax[j]    = 0.5 * (fMax[j] - fMin[j]);
    }

    //fprintf(err, "crop(): fMin=(%f %f %f) fMax=(%f %f %f) fCenter=(%f %f %f)\n",
    //    fMin[0], fMin[1], fMin[2],
    //    fMax[0], fMax[1], fMax[2],
    //    node->bnd.fCenter[0],
    //    node->bnd.fCenter[1],
    //    node->bnd.fCenter[2]);

}
#undef EXTEND_BOUNDS

static void split(pTree *node, int *d, FLOAT *fSplit)
{
    if (node->bnd.fMax[0] < node->bnd.fMax[1])
        *d = node->bnd.fMax[1] < node->bnd.fMax[2] ? SPLIT_Z_DIM : SPLIT_Y_DIM;
    else
        *d = node->bnd.fMax[0] < node->bnd.fMax[2] ? SPLIT_Z_DIM : SPLIT_X_DIM;

    *fSplit = node->bnd.fCenter[*d];

    //fprintf(err, "split(): %i %f\n", *d, *fSplit);
}


#define NEW_NODE(parent, child, l, u) \
{\
    parent->child = alloc_node(); \
    parent->child->pLeft = parent->child->pRight = NULL; \
    parent->child->pParent = parent; \
    parent->child->iLower = l; \
    parent->child->iUpper = u; \
}

pTree **sortStack = NULL;     /* this is the stack */

int partitionCount;

static int sort(struct pqNode **ps, pTree *node, int M)
{
    int sp=0;
    pTree *pLeft=NULL, *pRight=NULL;
    int nBucket = 0;

    /*========================================================================
     * Allocate stack
     *======================================================================*/
    int ns = (int)MMAX(1, floor(log(((double)(NPARTICLES+1))/(M+1))/log(2.0)));
    sp = 0;
    sortStack = (pTree **)realloc(sortStack, ns * sizeof(pTree *));
    assert(sortStack != NULL);

    if (node->iUpper - node->iLower + 1 <= M) 
        return 0;

    while (1)
    {
        int i, j;
        struct pqNode *t;
        int nl, nr;
        int d;
        FLOAT fSplit;

        assert(node != NULL);
        split(node, &d, &fSplit);

        i = node->iLower;
        j = node->iUpper;

        partitionCount++;
        PARTITION(ps, t, ->r[d], i,j, < fSplit,> fSplit);
        nl = i - node->iLower;
        nr = node->iUpper - i + 1;

        //fprintf(err, "nl=%i nr=%i\n", nl, nr);

        /*========================================================================
         * If both sides of the partition are not empty then create two new nodes
         * and crop the bounding boxes. Otherwise, undo the split operation.
         *======================================================================*/
        if (nl > 0 || nr > 0) {  // jpc changed this from && to || so that a split is always created.
            NEW_NODE(node, pLeft, node->iLower, i - 1);
            NEW_NODE(node, pRight, i, node->iUpper);
            nodeCount += 2;
#if 0
            crop(ps, node->pLeft, node->pLeft->iLower,  node->pLeft->iUpper);
            crop(ps, node->pRight, node->pRight->iLower, node->pRight->iUpper);
#else

            _3x_(i) node->pLeft->bnd.fCenter[i]  = node->bnd.fCenter[i];
            _3x_(i) node->pLeft->bnd.fMax[i]     = node->bnd.fMax[i];

            _3x_(i) node->pRight->bnd.fCenter[i] = node->bnd.fCenter[i];
            _3x_(i) node->pRight->bnd.fMax[i]    = node->bnd.fMax[i];

            node->pLeft->bnd.fMax[d]            *= 0.5;
            node->pLeft->bnd.fCenter[d]         -= node->pLeft->bnd.fMax[d];

            node->pRight->bnd.fMax[d]           *= 0.5;
            node->pRight->bnd.fCenter[d]        += node->pRight->bnd.fMax[d];
#endif

            pLeft = node->pLeft;
            pRight = node->pRight;
        }
        else
        {
            node->bnd.fMax[d] *= 0.5;
            if (nl > 0) {
                node->bnd.fCenter[d] -= node->bnd.fMax[d];
                pLeft = node;
            }
            else {
                node->bnd.fCenter[d] += node->bnd.fMax[d];
                pRight = node;
            }
        }

        /*========================================================================
         * Now figure out which subfile to process next
         *======================================================================*/
        if (nl > M && nr > M) 
        {
            if (nr > nl) sortStack[sp++] = pRight, node = pLeft;
            else sortStack[sp++] = pLeft, node = pRight;
        }
        else 
        {
            if (nl > M) node = pLeft;
            //else if (nl > 0) pLeft->iLower = 0;

            if (nr > M) node = pRight;
            //else if (nr > 0) pRight->iLower = 0;
        }

        if (nl <= M && nr <= M) {
            if (sp) node = sortStack[--sp];      /* pop tn */
            else break;
        }

    }

    return 0;
}

static void populate(struct pqNode **ps, pTree * const node)
{
    FLOAT m, fMass, x, y, z;

    if (node == NULL) return;

    populate(ps, node->pLeft);
    populate(ps, node->pRight);

    m = fMass = x = y = z = 0;

    if (node->pLeft == NULL  &&  node->pRight == NULL)
    {
        /*====================================================================
         * Calculate total mass and CoM for a leaf node
         *==================================================================*/
        int pj = node->iLower;
        //ps[pj].iBucket = iNode;
        for (pj = node->iLower; pj <= node->iUpper; ++pj) 
        {
            //p[pj].iBucket = iNode;
            m = ps[pj]->fMass;
            fMass += m;
            x += m * ps[pj]->r[0];
            y += m * ps[pj]->r[1];
            z += m * ps[pj]->r[2];
        }
        m = 1/fMass;
        node->r[0]  = m * x;
        node->r[1]  = m * y;
        node->r[2]  = m * z;
        node->fMass = fMass;
    }
    else
    {
        /*====================================================================
         * Calculate total mass and CoM for a node
         *==================================================================*/
        pTree *pLeft = node->pLeft;
        pTree *pRight = node->pRight;

        if (pLeft != NULL)
        {
            fMass = pLeft->fMass;
            x = pLeft->fMass * pLeft->r[0];
            y = pLeft->fMass * pLeft->r[1];
            z = pLeft->fMass * pLeft->r[2];
        }

        if (pRight != NULL)
        {
            fMass += pRight->fMass;
            x += pRight->fMass * pRight->r[0];
            y += pRight->fMass * pRight->r[1];
            z += pRight->fMass * pRight->r[2];
        }

        m = 1/fMass;
        node->r[0]  = m * x;
        node->r[1]  = m * y;
        node->r[2]  = m * z;
        node->fMass = fMass;
    }
}

pTree *pintBuildTree(struct pqNode **ps, pTree *head)
{
    nodePoolUsed = 0;
    memset(node_pool, 0, nodePoolSize * sizeof(pTree));

    head = NULL;
    if (head == NULL)
    {
        head = alloc_node();
        assert(head != NULL);
        head->iLower = 0;
        head->iUpper = NPARTICLES-1;
    }

    comparisonCount = 0;
    nodeCount = 1;

    crop(ps, head, 0, NPARTICLES-1);

    //fprintf(err, "(%e %e %e) (%e %e %e)\n", head->bnd.fMax[0], head->bnd.fMax[1], head->bnd.fMax[2], head->bnd.fCenter[0], head->bnd.fCenter[1], head->bnd.fCenter[2]);

    sort(ps, head, MAX_PIC);
    //populate(ps, head);

    _DA_ fprintf(err, "nodePoolUsed = %i\n", nodePoolUsed);
    _DA_ fprintf(err, "head->fMass = %f\n", head->fMass);
    _DA_ fprintf(err, "head->r = %f %f %f\n", head->r[0], head->r[1], head->r[2]);


    return head;
}

#else /* DYNAMIC_TREE */

static void split2(pTree *node, int *d, FLOAT *fSplit)
{
    if (node->cell.fMax[0] < node->cell.fMax[1])
        *d = node->cell.fMax[1] < node->cell.fMax[2] ? SPLIT_Z_DIM : SPLIT_Y_DIM;
    else
        *d = node->cell.fMax[0] < node->cell.fMax[2] ? SPLIT_Z_DIM : SPLIT_X_DIM;

    *fSplit = node->cell.r[*d];

    //fprintf(err, "split(): %i %f\n", *d, *fSplit);
}

double inbounds_time;
int inbounds_ncalls;

int inbounds(struct pqNode *p, pTree *node)
{
    //_DA_ fprintf(err, "inbounds(): %f %f %f\n", p->r[0], p->r[1], p->r[2]);
    //_DA_ fprintf(err, "            %f %f %f\n", cell->r[0], cell->r[1], cell->r[2]);
    //_DA_ fprintf(err, "            %f %f %f\n", cell->fMax[0], cell->fMax[1], cell->fMax[2]);

    //inbounds_ncalls++;

    return (FABS(p->r[0] - node->cell.r[0]) <= node->cell.fMax[0] + 1e-5) 
        && (FABS(p->r[1] - node->cell.r[1]) <= node->cell.fMax[1] + 1e-5) 
        && (FABS(p->r[2] - node->cell.r[2]) <= node->cell.fMax[2] + 1e-5);
}

#if 1
static void assign_split(TREE_TYPE node, int *d0, FLOAT *fSplit0)
{
    int d;

    if (NODE(node)->cell.fMax[0] < NODE(node)->cell.fMax[1])
        d = NODE(node)->cell.fMax[1] < NODE(node)->cell.fMax[2] ? SPLIT_Z_DIM : SPLIT_Y_DIM;
    else
        d = NODE(node)->cell.fMax[0] < NODE(node)->cell.fMax[2] ? SPLIT_Z_DIM : SPLIT_X_DIM;

    assert(d != 3);
    NODE(node)->cell.fSplit = NODE(node)->cell.r[d];
    NODE(node)->cell.d = d;
    //SPLIT_FSPLIT(NODE(node)->cell.split) = NODE(node)->cell.r[d];
    //SET_SPLIT_D(NODE(node), d);

    /**/ //assert(SPLIT_D(NODE(node)->cell.split) == d);

    *d0 = d;
    *fSplit0 = NODE(node)->cell.r[d];
}
#endif

struct pqNode *join_lists(struct pqNode *a, struct pqNode *b)
{
    struct pqNode *prev, *t;

    if (a == NULL) return b;
    if (b == NULL) return a;

    for (prev = t = a; t != NULL; prev = t, t = t->neighbor) {}

    prev->neighbor = b;

    return a;
}

#if 1
#if 0
void split_list(pTree *head, struct pqNode **left, struct pqNode **right, int *lcount, int *rcount)
{
    struct pqNode *p;
    struct pqNode *leftList  = NULL, 
                  *rightList = NULL;

    int d = head->d;
    FLOAT fSplit = head->fSplit;

    int nl=0, nr=0;

    for (p = head->cell.list; p != NULL; p = p->neighbor)
    {
        comparisonCount++;

        if (p->r[d] < fSplit)
        {
            p->node = head->pLeft;

            if (leftList == NULL) leftList = *left = p;
            else                  leftList = leftList->neighbor     = p;

            /**/ assert(p->node->cell.list != NULL);

            //head->pLeft->cell.count++;
            nl++;
        }
        else
        {
            p->node = head->pRight;

            if (rightList == NULL) rightList = *right = p;
            else                   rightList = rightList->neighbor     = p;

            /**/ assert(p->node->cell.list != NULL);

            //head->pRight->cell.count++;
            nr++;
        }
    }

    *lcount = nl;
    *rcount = nr;

    if (leftList  != NULL) leftList->neighbor  = NULL;
    if (rightList != NULL) rightList->neighbor = NULL;
}
#else
void split_list(long head, struct pqNode **left, struct pqNode **right, int *lcount, int *rcount)
{
    struct pqNode *p;
    struct pqNode *leftList  = NULL, 
                  *rightList = NULL;

    int d = NODE(head)->cell.d; //SPLIT_D(head->cell.split); //head->cell.split.d;
    FLOAT fSplit = NODE(head)->cell.fSplit;

    int nl=0, nr=0;

#if 1
    if (NODE(head)->cell.left_count == 0)
    {
        *right  = NODE(head)->cell.list;
        *rcount = NODE(head)->cell.count;
        return;
    }

    if (NODE(head)->cell.left_count == NODE(head)->cell.count)
    {
        *left   = NODE(head)->cell.list;
        *lcount = NODE(head)->cell.count;
        return;
    }
#endif

    for (p = NODE(head)->cell.list; p != NULL; p = p->neighbor)
    {
        comparisonCount++;

        if (p->r[d] < fSplit)
        {
            p->node = NODE(NODE(head)->pLeft);

            if (leftList == NULL) leftList = *left = p;
            else                  leftList = leftList->neighbor     = p;

            /**/ //assert(p->node->cell.list != NULL);

            //NODE(head)->pLeft->cell.count++;
            nl++;
        }
        else
        {
            p->node = NODE(NODE(head)->pRight);

            if (rightList == NULL) rightList = *right = p;
            else                   rightList = rightList->neighbor     = p;

            /**/ //assert(p->node->cell.list != NULL);

            //NODE(head)->pRight->cell.count++;
            nr++;
        }
    }

    *lcount = nl;
    *rcount = nr;

    if (leftList  != NULL) leftList->neighbor  = NULL;
    if (rightList != NULL) rightList->neighbor = NULL;
}
#endif

#else
void split_list(pTree *head)
{
    struct pqNode *p;
    struct pqNode *leftList  = NULL, 
                  *rightList = NULL;

    int d = head->d;
    FLOAT fSplit = head->fSplit;

    struct pqNode *sp, *ep, *prev;

    sp = head->cell.list;
    assert(head->cell.count <= MAX_PIC);
    assert(head->cell.count > 0);

    _DA_ fprintf(err, "head->cell.count=%i\n", head->cell.count);

    if (sp == NULL) return;

    int count = 0;
    for (prev = ep = sp; ep != NULL; ep = ep->neighbor) 
    {
        comparisonCount++;
        if (ep->r[d] >= fSplit) break;
        count++;
        prev = ep;
    }

    if (count == 0)
    {
        head->pRight->cell.list = sp;
        head->pRight->cell.count = head->cell.count;
        return;
    }
    else if (count == head->cell.count)
    {
        head->pLeft->cell.list = sp;
        head->pLeft->cell.count = head->cell.count;
        return;
    }

    while (sp != NULL)
    {
    comparisonCount--;
        count = 0;
        for (prev = ep = sp; ep != NULL; ep = ep->neighbor)
        {
            comparisonCount++;
            _DA_ fprintf(err, "ep->r[d]=%f fSplit=%f count=%i\n", ep->r[d], fSplit, count);
            if (ep->r[d] < fSplit) break;
            count++;
            prev = ep;
        }

        if (count)
        {
            if (rightList == NULL) rightList = head->pRight->cell.list = sp;
            else                   rightList->neighbor = sp;
            head->pRight->cell.count += count;

            sp = prev->neighbor;
            prev->neighbor = NULL;
            rightList = prev;
        }
    comparisonCount--;

        count = 0;
        for (prev = ep = sp; ep != NULL; ep = ep->neighbor) 
        {
            comparisonCount++;
            _DA_ fprintf(err, "ep->r[d]=%f fSplit=%f count=%i\n", ep->r[d], fSplit, count);
            if (ep->r[d] >= fSplit) break;
            count++;
            prev = ep;
        }

        if (count)
        {
            if (leftList == NULL) leftList = head->pLeft->cell.list = sp;
            else                  leftList->neighbor = sp;
            head->pLeft->cell.count += count;

            sp = prev->neighbor;
            prev->neighbor = NULL;
            leftList = prev;
        }
    }

}
#endif

#if 1
//if (head->T != NULL) { head->T->cell.list = NULL; head->T->cell.count= head->T->cell.left_count = 0;  } \
//else \

    //SET_SPLIT_D(head->T, 3); \

#define CREATE_BRANCH(head, d, T, OP) \
{\
    NODE(head)-> T  = alloc_node(); NODE(NODE(head)-> T) ->pParent  = head;      \
    NODE(NODE(head)-> T) -> cell.d = 4; \
    _3x_(i)                                                    \
    {                                                          \
        NODE(NODE(head)-> T) ->cell.r[i]      = NODE(head)->cell.r[i];           \
        NODE(NODE(head)-> T) ->cell.fMax[i]   = NODE(head)->cell.fMax[i];        \
    }                                                          \
    NODE(NODE(head)-> T) ->cell.fMax[d]      *= 0.5;                       \
    NODE(NODE(head)-> T) ->cell.r[d]         = NODE(NODE(head)-> T) ->cell.r[d] OP NODE(NODE(head)-> T) ->cell.fMax[d];   \
    nodeCount++; \
}
    //assign_split(head-> T );                                   \

#else
#define CREATE_BRANCH(head, d, T, OP) \
{\
    if (head->T != NULL) { /* head->T->cell.list = NULL; head->T->cell.count= head->T->cell.left_count = 0; */ } \
    else { \
    head-> T  = alloc_node(); head-> T ->pParent  = head;      \
    head->has_children = 1; \
    switch (d) { \
        case 0:  \
            head-> T ->cell.fMax[1] = head->cell.fMax[1]; \
            head-> T ->cell.fMax[2] = head->cell.fMax[2]; \
            head-> T ->cell.fMax[0] = 0.5 * head->cell.fMax[0];          \
            head-> T ->cell.r[0]    = head->cell.r[0] OP head-> T ->cell.fMax[0]; \
            head-> T ->cell.r[1]    = head->cell.r[1]; \
            head-> T ->cell.r[2]    = head->cell.r[2]; \
            break; \
        case 1:  \
            head-> T ->cell.fMax[0] = head->cell.fMax[0]; \
            head-> T ->cell.fMax[1] = 0.5 * head->cell.fMax[1];          \
            head-> T ->cell.fMax[2] = head->cell.fMax[2]; \
            head-> T ->cell.r[0]    = head->cell.r[0]; \
            head-> T ->cell.r[1]    = head->cell.r[1] OP head-> T ->cell.fMax[1]; \
            head-> T ->cell.r[2]    = head->cell.r[2]; \
            break; \
        case 2:  \
            head-> T ->cell.fMax[0] = head->cell.fMax[0]; \
            head-> T ->cell.fMax[1] = head->cell.fMax[1]; \
            head-> T ->cell.fMax[2] = 0.5 * head->cell.fMax[2];          \
            head-> T ->cell.r[0]    = head->cell.r[0]; \
            head-> T ->cell.r[1]    = head->cell.r[1]; \
            head-> T ->cell.r[2]    = head->cell.r[2] OP head-> T ->cell.fMax[2]; \
            break; \
    } \
    assign_split(head-> T);                                   \
    nodeCount++; \
    }\
}
    
#endif

#if 0
inline int split_cell(pTree *head)
{
    int i;
    int   d      = head->d;
    FLOAT fSplit = head->fSplit;

    CREATE_BRANCH(pLeft,  -);
    CREATE_BRANCH(pRight, +);

    /**/ assert(head->pLeft->cell.count == 0);
    /**/ assert(head->pRight->cell.count == 0);

    split_list(head);

    //_DA_ fprintf(err, "%e %e %e\n", head->cell.fMass, head->pLeft->cell.fMass, head->pRight->cell.fMass);
    //_DA_ fprintf(err, "%i %e\n", head->d, head->fSplit);
    _DA_ fprintf(err, "Left count = %i  Right count = %i\n", head->pLeft->cell.count, head->pRight->cell.count);
    // /**/ assert(head->cell.fMass == (head->pLeft->cell.fMass + head->pRight->cell.fMass));

    head->cell.list  = NULL;
    head->cell.count = 0;

    return d;
}
#else

void create_left_branch(long head, const int d)
{
        int i;
        CREATE_BRANCH(head, d, pLeft,  -);
}

void create_right_branch(long head, const int d)
{
        int i;
        CREATE_BRANCH(head, d, pRight,  +);
}

long split_cell(const struct pqNode *p, long head)
{
    int i;
    int lcount=0, rcount=0;
    int d;
    FLOAT fSplit;
    union split_u split;

    for (;;)
    {
        struct pqNode *left=NULL, *right=NULL;
        lcount = rcount = 0;

        assign_split(head, &d, &fSplit);

#if 0
        { 
            split  = head->cell.split; 
            d      = SPLIT_D(split); 
            fSplit = SPLIT_FSPLIT(split); 
        }
#endif

        /**/ assert(d != 4);

        split_list(head, &left, &right, &lcount, &rcount);

        _DA_ fprintf(err, "Left count = %i  Right count = %i\n", lcount, rcount);
        //_DA_ fprintf(err, "d=%i fSplit=%f\n", SPLIT_D(head->cell.split), head->cell.split.fSplit);

        NODE(head)->cell.list  = NULL;
        NODE(head)->cell.count = 0;
        NODE(head)->cell.left_count = 0;

        create_left_branch(head, d);
        NODE(NODE(head)->pLeft)->cell.list = left;
        NODE(NODE(head)->pLeft)->cell.count = lcount;

        create_right_branch(head, d);
        NODE(NODE(head)->pRight)->cell.list  = right;
        NODE(NODE(head)->pRight)->cell.count = rcount;

        comparisonCount++;
        if (p->r[d] < fSplit)
        {
            head = NODE(head)->pLeft;
            if (rcount != 0) break;
        }
        else
        {
            head = NODE(head)->pRight;
            if (lcount != 0) break;
        }

    }

    return head;
}

#endif

#define FLOAT2UINTCAST(f) (*((unsigned int *)(&f)))

//pTree *find_leaf(struct pqNode *p, register pTree *head)
TREE_TYPE find_leaf(FLOAT *r0, register long head)
{
    register FLOAT fSplit;
    register int d;
    //register pTree *hL;
    register TREE_TYPE hL;

    FLOAT r[3];

    r[0] = r0[0];
    r[1] = r0[1];
    r[2] = r0[2];

    for (;;) 
    {
        //union split_u split = head->cell.split; 
        d = NODE(head)->cell.d; //SPLIT_D(split); 

        if (d == 3) return head;

        fSplit = NODE(head)->cell.fSplit; //SPLIT_FSPLIT(split); 

        comparisonCount++;

        _DA_ fprintf(err, "d=%i fSplit=%f\n", d, fSplit);
 
        hL = NODE(head)->pLeft;

        head = NODE(head)->pRight;
        if (r[d] < fSplit) head = hL;

        assert(head != NULL);
    }

    return head;
}

TREE_TYPE tree_build(struct pqNode **ps, long old_head)
{
    /**/ assert(head != NULL); 
    _DA_ fprintf(err, "Finding leaf...\n");
    register long head;
    int i;

    double t=0;
    double start, end;

    for (i=0; i < NPARTICLES; i++)
    {
        register FLOAT fSplit;
        register int d;
        register long hL;
        register struct pqNode *p = ps[i];

        head = old_head;

        for (;;) 
        {
            d = NODE(head)->cell.d;

            if (d & 0x4) break;

            fSplit = NODE(head)->cell.fSplit;

            comparisonCount++;

            hL = NODE(head)->pLeft;

            head = NODE(head)->pRight;
            if (p->r[d] < fSplit) head = hL;
        }

        //start = CPUTIME;
        if (NODE_FULL(NODE(head))) 
        {
            _DA_ fprintf(err, "Splitting cell...\n");
            head = split_cell(p, head);
        }
        //end = CPUTIME;
        //t += (end - start);

        d      = NODE(head)->cell.d; //SPLIT_D(split);
        fSplit = NODE(head)->cell.fSplit; //split.fSplit;

        if (p->r[d] < fSplit)
            NODE(head)->cell.left_count++;

        p->node     = NODE(head);
        p->neighbor = NODE(head)->cell.list;

        NODE(head)->cell.list = p;
        NODE(head)->cell.count++;
    }

    //fprintf(err, "t=%2.5e\n", t);
    return old_head;
}

#if 0
TREE_TYPE tree_insert(struct pqNode *p, TREE_TYPE head)
{
    /**/ assert(p != NULL);
    /**/ assert(head != NULL);

    //double start = CPUTIME;
    //bottom_up_ncalls++;

    int i;

    while (1)
    {
        pTree *old_head = head;
            
        //if (inbounds(p, head)) // by construction of the head, we don't need this for building the tree
        {
            /**/ assert(head != NULL); 
            //head = find_leaf(p, head);
            _DA_ fprintf(err, "Finding leaf...\n");
            //head = find_leaf(p->r, head);


    register FLOAT fSplit;
    register int d;
    register pTree *hL;
    for (;;) 
    {
        d = NODE(head)->cell.d;

        if (d == 3) break;

        fSplit = head->cell.fSplit;

        comparisonCount++;

        hL = NODE(head)->pLeft;

        head = NODE(head)->pRight;
        if (p->r[d] < fSplit) head = hL;
    }


            if (NODE_FULL(head)) 
            {
                _DA_ fprintf(err, "Splitting cell...\n");
                head = split_cell(p, head);
            }

            assert(head != NULL);

            //union split_u split = head->cell.split;
            d      = head->cell.d; //SPLIT_D(split);
            fSplit = head->cell.fSplit; //split.fSplit;

            if (p->r[d] < fSplit)
                head->cell.left_count++;

            p->node     = head;
            p->neighbor = head->cell.list;

            head->cell.list = p;
            head->cell.count++;
            assert(head->cell.count <= MAX_PIC);

            return old_head;
        }

        comparisonCount++;

        if (head->pParent == NULL) break;
        head = head->pParent;
    }

    _DA_ fprintf(err, "Creating new parent/child [id=%i]\n", ID(p));
    int dir = 1;
    int d;
    int best_dir=1, best_dim = SPLIT_X_DIM;
    float best_dist=1e38;
    float best_edge = 0;

    /*================================================================
     * Find the face that is closest to the particle.
     * (Could swap the loops and change dir = -1 to dir *= -1)
     *==============================================================*/
    for (i=0; i < 2; i++) /* two directions */
    {
        for (d=0; d < 3; d++) /* three dimensions */
        {
            float edge = head->cell.r[d] + head->cell.fMax[d] * dir;
            float dist = p->r[d] - edge;

            if (FABS(dist) < FABS(best_dist))
            {
                best_dist = dist;
                best_dir  = dir;
                best_dim  = d;
                best_edge = edge;
            }
        }

        dir = -1;
    }

    //double end = CPUTIME;
    //_DA_ bottom_up_time += (end - start);

    _DA_ fprintf(err, "    best_dim=%i best_dir=%i\n", best_dim, best_dir);

    pTree *newHead  = alloc_node();
    pTree *newChild = alloc_node();

    _3x_(i)
    {
        newHead->cell.r[i]     = head->cell.r[i];
        newChild->cell.r[i]    = head->cell.r[i];

        newHead->cell.fMax[i]  = head->cell.fMax[i];
        newChild->cell.fMax[i] = head->cell.fMax[i];
    }

    newChild->pParent              = newHead;
    newChild->pLeft                = NULL;
    newChild->pRight               = NULL;
    //newChild->cell.fMass           = 0;
    newChild->cell.r[best_dim]    += 2 * head->cell.fMax[best_dim] * best_dir;
    newChild->cell.count           = 0;
    newChild->cell.list            = NULL;

    newHead->pParent               = NULL;
    //newHead->cell.fMass            = head->cell.fMass;
    newHead->cell.r[best_dim]     += head->cell.fMax[best_dim] * best_dir;
    newHead->cell.fMax[best_dim]  *= 2;
    newHead->cell.count            = 0;
    newHead->cell.list             = NULL;

    //assign_split(newHead);
    //assign_split(newChild);

    _DA_ fprintf(err, "(%e) %e %e %e\n", 
        head->cell.fMax[best_dim], head->cell.r[best_dim], newChild->cell.r[best_dim], newHead->cell.r[best_dim]);

    comparisonCount++;

    if (best_dir == -1)
    {
        newHead->pLeft  = newChild;
        newHead->pRight = head;
    }
    else
    {
        newHead->pLeft  = head;
        newHead->pRight = newChild;
    }

    /**/ //assert(HAS_CHILDREN(newHead));

    //_DA_ fprintf(err, "Particle still out of bounds [id=%i]\n", ID(p));
    return tree_insert(p, newHead);
}
#endif

#if 0
pTree *tree_remove(struct pqNode *p, pTree *head)
{
    _DA_ fprintf(err, "remove() [id=%i]\n", ID(p));

    /**/ assert (p->node != NULL);
    /**/ assert (p->node->cell.list != NULL);

    pTree *node = p->node;
    struct pqNode *list = node->cell.list;
    struct pqNode *item, *temp;

    /*========================================================================
     * Find the particle in the list. item will point to the correct place on 
     * exit. Notice that there is no check for loop termination other than the 
     * assert. It is a serious error if the particle is not in the list.
     *======================================================================*/
    for (item = temp = list; item != p; temp = item, item = item->neighbor)
        assert(item != NULL);

    if (item == list) node->cell.list = list->neighbor; /* first item in list */
    else              temp->neighbor  = item->neighbor;

    node->cell.count--;                 /**/ assert(node->cell.count >= 0);

#if 0
    for (; node != NULL; node = node->pParent)
    {
        node->cell.fMass -= MASS(p);    /**/ assert(node->cell.fMass >= -1e-3);
    }
#endif

    node = p->node;

    pTree *parent = node->pParent;

    /*========================================================================
     * If we are not at the root node then examine the sibling to see if we
     * can merge the children it back into the parent.
     *======================================================================*/
    if (parent != NULL)
    {
        /**/ assert(parent->cell.count == 0);
        /**/ assert(HAS_CHILDREN(parent));

        pTree *sibling = (parent->pLeft == node) ? parent->pRight : parent->pLeft;

        /* node doesn't have children because p is in node */ assert(NO_CHILDREN(node));   

        if (NO_CHILDREN(sibling)) 
        {
            if (node->cell.count + sibling->cell.count <= MAX_PIC)
            {
                _DA_ fprintf(err, "Merging children [id=%i]\n", ID(p));

                //parent->cell.fMass = sibling->cell.fMass + node->cell.fMass;
                parent->cell.count = sibling->cell.count + node->cell.count;
                parent->cell.list  = join_lists(sibling->cell.list, node->cell.list);

                int c=0;
                for (temp = parent->cell.list; temp != NULL; temp = temp->neighbor)
                {
                    temp->node = parent;
                    c++;
                }

                /**/ assert(c == node->cell.count + sibling->cell.count);

                parent->pLeft      = NULL;
                parent->pRight     = NULL;

                if (head == node) head = parent;

                free_node(sibling); sibling = NULL;
                free_node(node);    node = NULL;
            }
        }
    }

    p->node = NULL;

    free(item); 

    return head;
}
#endif

#define EXTEND_BOUNDS(R, MIN, MAX) \
{      if (R < MIN) MIN = R; \
  else if (R > MAX) MAX = R; }
static void crop(struct pqNode **ps, TREE_TYPE node, int front, int back)
{
    int j;
    FLOAT fMin[3],fMax[3];

    //fprintf(err, "crop(): %i %i\n", front, back);

    _3x_(j)
    {
        fMin[j] = ps[front]->r[j];
        fMax[j] = ps[front]->r[j];
    }

    while (++front <= back)
    {
        EXTEND_BOUNDS(ps[front]->r[0], fMin[0], fMax[0]);
        EXTEND_BOUNDS(ps[front]->r[1], fMin[1], fMax[1]);
        EXTEND_BOUNDS(ps[front]->r[2], fMin[2], fMax[2]);
    }

    _3x_(j)
    {
        NODE(node)->cell.r   [j] = 0.5 * (fMax[j] + fMin[j]);
        NODE(node)->cell.fMax[j] = 0.5 * (fMax[j] - fMin[j]);
        assert(NODE(node)->cell.fMax[j] > 0);
    }

    //fprintf(err, "crop(): fMin=(%f %f %f) fMax=(%f %f %f) fCenter=(%f %f %f)\n",
    //    fMin[0], fMin[1], fMin[2],
    //    fMax[0], fMax[1], fMax[2],
    //    node->bnd.fCenter[0],
    //    node->bnd.fCenter[1],
    //    node->bnd.fCenter[2]);

}
#undef EXTEND_BOUNDS

pTree *pintBuildTree(struct pqNode **ps, pTree *head)
{
    int i=0;
    nodePoolUsed = 0;
    itemPoolUsed=0;
    TREE_TYPE h;

    memset(node_pool, 0, nodePoolSize * sizeof(pTree));

    if (head == NULL)
    {
        //_DA_ fprintf(err, "Creating head: cellsize=%f\n", cellSize);
        h = alloc_node();
        NODE(h)->cell.d = 4;
        //SET_SPLIT_D(head, 3);
    }
    else
    {
        //h = NODE(head);
    }

    crop(ps, h, 0, NPARTICLES-1);
    //assign_split(head);

    comparisonCount = 0;
    nodeCount = 1;

    h = tree_build(ps, h);
#if 0
    for (i=0; i < NPARTICLES; i++)
    {
        head = tree_build(ps[i], head);
        //head = tree_insert(ps[i], head);
    }
#endif

    return NODE(h);
}

#endif /* DYNAMIC_TREE */

#if 0
int main(int argc, char **argv)
{
    int i, t, j;
    int a[10];
    for (i=0; i < 10; i++)
        a[i] = 10-i;

    i=0;
    j=9;
    PARTITION(a, t, /*EMPTY*/, i,j, % 2, %2 == 0);

    for (i=0; i < 10; i++)
        printf("%i\n", a[i]);


    return 0;
}
#endif

