#ifndef PARTICLE_INCLUDED
#define PARTICLE_INCLUDED

typedef float real;

struct particle {
    real r[3];
    real v[3];
    real a[3];

    int id;
    unsigned short iType;
    unsigned char pqDist;       /* can be small, N < 2^256 should be safe */
    real mass, soft;
    char pad[12];
};

#define MASS(p) ((p)->mass)
#define SOFT(p) ((p)->soft)
#define ID(p) ((p)->id)
#define VEL(p) ((p)->v)
#define RADIUS(p) ((p)->r)
#define POS(p) (RADIUS(p))

#define PQ_INIT(p) {\
    (p)->pqDist = 1;\
    }

#define PQ_MERGE(p,q) {\
    struct pqNode *PQ_t=NULL,*PQ_r = NULL; unsigned char PQ_d=0,PQ_dt=0; \
    while ((p) && (q)) { \
        if ((p)->timeNext < (q)->timeNext) { \
            PQ_t = (p)->pqRight; \
            (p)->pqRight = PQ_r; \
            PQ_r = (p); \
            (p) = PQ_t; \
            } \
        else { \
            PQ_t = (q)->pqRight; \
            (q)->pqRight = PQ_r; \
            PQ_r = (q); \
            (q) = PQ_t; \
            } \
        } \
    if (!(p)) (p) = (q); \
    if ((p)) PQ_d = (p)->pqDist; /* only a merge of two null trees needs this if */ \
    while (PQ_r) { \
        (q) = PQ_r->pqRight; \
        PQ_t = PQ_r->pqLeft; \
        if (PQ_t) PQ_dt = PQ_t->pqDist; \
        else PQ_dt = 0; \
        if (PQ_dt < PQ_d) { \
            PQ_d = PQ_dt + 1; \
            PQ_r->pqRight = PQ_t; \
            PQ_r->pqLeft = (p); \
            } \
        else { \
            ++PQ_d; \
            PQ_r->pqRight = (p); \
            } \
        PQ_r->pqDist = PQ_d; \
        (p) = PQ_r; \
        PQ_r = (q); \
        } \
    }

#define PQ_MIN(head) (head)

#define PQ_REMOVE_MIN(head, p) { \
    struct pqNode *PQ_l = (head)->pqLeft; \
    struct pqNode *PQ_r = (head)->pqRight; \
    struct pqNode *t; \
    (p) = head; \
    (p)->pqLeft = NULL; (p)->pqRight = NULL; (p)->pqDist = 1;\
    head = PQ_l; \
    t = PQ_r; \
    PQ_MERGE(head, t); \
}

#endif
