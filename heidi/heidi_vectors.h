#ifndef HEIDI_VECTORS
#define HEIDI_VECTORS

#define NVECTORELEMENTS 3

#define NOP ((void)(0))
#define _ [__]
#define Vector(type, var) type var[NVECTORELEMENTS]
#define V(op) VEC(NOP, op, NOP)
#define VECELEM(vec, elem) (vec [ elem ])

#if NVECTORELEMENTS == 3

#define NORM(v) SQRT((v)[0]*(v)[0] + (v)[1]*(v)[1] + (v)[2]*(v)[2])
#define UNIT(v, o) do { __typeof__ ((v)[0]) norm = NORM(v); V( (o)_ = v _ / norm ); } while(0)

#define VEC(preamble, op, postamble) \
do {\
preamble;   \
{ const unsigned char __ = 0; do { op; } while(0); } \
{ const unsigned char __ = 1; do { op; } while(0); } \
{ const unsigned char __ = 2; do { op; } while(0); } \
postamble;  \
} while(0)

#define VECIF(op) \
{\
int vecif_res = 1; \
{ const unsigned char __ = 0; vecif_res = vecif_res && (op); } \
{ const unsigned char __ = 1; vecif_res = vecif_res && (op); } \
{ const unsigned char __ = 2; vecif_res = vecif_res && (op); } \
if (vecif_res) {

#define ELSEVECIF } else {
#define ENDVECIF }} 

#define CROSS(u, v, q) \
do {\
    (q)[0] = (u)[1]*(v)[2] - (u)[2]*(v)[1]; \
    (q)[1] = (u)[0]*(v)[2] - (u)[2]*(v)[0]; \
    (q)[2] = (u)[0]*(v)[1] - (u)[1]*(v)[0]; \
} while(0)

#define DOT(u, v) ((u)[0] * (v)[0] + (u)[1] * (v)[1] + (u)[2] * (v)[2])

#endif

#endif
