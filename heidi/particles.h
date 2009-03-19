#ifndef PARTICLES_H
#define PARTICLES_H

#include "heidi_types.h"


#if 0
enum ParticleType
{
    ELECTRON=0,
    PROTON,
    NEUTRON,
    PHOTON,
    NPARTICLETYPES,
};
#endif

typedef struct 
{
    Vector(Pos, x);
    Vector(Vel, v);
    Vector(Acc, a);

    Integer type;
    Time dt;

} Particle;

#endif

