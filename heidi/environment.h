#ifndef HEIDI_ENVIRONMENT_H
#define HEIDI_ENVIRONMENT_H

#include "heidi_types.h"
#include "io.h"

typedef struct
{
    Real    G;          /* Gravitational Constant  */
    Real    kC;         /* Coulomb's Constant      */
    Real    kB;         /* Boltzmann constant      */
    Speed   c;          /* Speed of light          */
    Real    fsc;        /* fine structure constant */
    Real    hbar;       /* h / 2pi                 */

    Mass    Me;         /* Electron mass */
    Mass    Mp;         /* Proton mass   */
    Mass    Mn;         /* Neutron mass  */

    Mass    Qe;         /* Electron charge */
    Mass    Qp;         /* Proton charge   */
    Mass    Qn;         /* Neutron charge  */
  
    Real    h;          /* Plank's constant */

} Constants;

typedef struct
{
    Integer N;
    Particle *p;
} ParticleList;

typedef struct
{
    Mass m;
    Charge q;
    Integer Z, N, E;
} ParticleTypeData;

typedef struct
{
    Real kin, pot;
} Energy;

typedef struct
{
    Constants K;
    Energy E;

    Time t_now;
    Time t_end;
    Integer t_step;

    IO io;

    ParticleTypeData *particle_type_data;
    ParticleList plist;

} Environment;

void init_env();

#endif
