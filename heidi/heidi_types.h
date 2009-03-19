#ifndef HEIDI_TYPES_H
#define HEIDI_TYPES_H

#include "heidi_vectors.h"

#define ABS(x)   fabs(x)
#define MIN(x,y) fmin(x,y)
#define MAX(x,y) fmax(x,y)
#define SQRT(x)  sqrt(x)
#define SWAP(x, y) do { __typeof__ (x) t; t = (x); (x) = (y); (y) = t; } while (0)

typedef double Real;
typedef int   Integer;

typedef Real Mass;
typedef Real Coord;
typedef Real Speed;
typedef Real Pos;
typedef Real Vel;
typedef Real Acc;
typedef Real Dist;
typedef double Time;
typedef Real Charge;

typedef Integer ID; 


typedef struct
{ 
    Integer len;
    ID  id[0]; 
} IDList;

typedef IDList ProtonList;
typedef IDList ElectronList;

#endif
