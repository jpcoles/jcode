#ifndef HEIDI_H
#define HEIDI_H

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "heidi_types.h"

#include "environment_macros.h"  
#include "particles.h"  
#include "physical_constants.h"
#include "environment.h"

#define PTR_MALLOC(T, s) ((T *)malloc((s) * sizeof(T)))
#define PTR_CALLOC(T, s) ((T *)calloc((s), sizeof(T)))

#define LIST_MALLOC(T1, s1, T2, s2) ((T1 *)malloc((s1) * (sizeof(T1) + (s2) * sizeof(T2))))

#define VL(x) if (env.io.verbosity >= (x))
#define STDOUT fprintf(env.io.out, 
#define STDERR fprintf(env.io.err, 

#define all_particles(i) ((i)=0; (i) < NP; (i)++)
#define all_particle_pairs(i,j) ((i)=0; (i) < NP-1; (i)++) for ((j)=(i)+1; (j) < NP; (j)++)

#endif

