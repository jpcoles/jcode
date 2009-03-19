#include <math.h>

#include "heidi.h"
#include "environment.h"

Environment env;

void default_constants()
{
    env.K.G    =   GRAVITATIONAL_CONSTANT;
    env.K.Me   =   ELECTRON_MASS;
    env.K.Mp   =   PROTON_MASS;
    env.K.Mn   =   NEUTRON_MASS;
    env.K.Qe   =   ELECTRON_CHARGE;
    env.K.Qp   =   PROTON_CHARGE;
    env.K.Qn   =   NEUTRON_CHARGE;

    env.K.h    =   PLANKS_CONSTANT_H;
    env.K.hbar =   PLANKS_CONSTANT_HBAR;

    env.K.c    =   SPEED_OF_LIGHT;

    env.K.fsc  =   FINE_STRUCTURE_CONSTANT;
    env.K.kC   =   COULOMB_CONSTANT;
    env.K.kB   =   BOLTZMANN_CONSTANT;
}

void default_io()
{
    env.io.in  = stdin;
    env.io.out = stdout;
    env.io.err = stderr;

    env.io.state_fp = NULL;
    env.io.stats_fp = NULL;

    env.io.state_file[0] = '\0';
    env.io.verbosity = 100;
}

#if 0
void default_particle_type_data()
{
    env.particle_type_data[ELECTRON].m = 4*_Mp;
    env.particle_type_data[ELECTRON].q = 2*_Qp;

    env.particle_type_data[PROTON].m   = 2*79*_Mp;
    env.particle_type_data[PROTON].q   = 79*_Qp;

    env.particle_type_data[NEUTRON].m  = _Mn;
    env.particle_type_data[NEUTRON].q  = _Qn;

    env.particle_type_data[PHOTON].q   = 0;
    env.particle_type_data[PHOTON].m   = 0;
}
#endif

void init_env()
{
    default_constants(); /* This should come first */
    default_io();
    //default_particle_type_data();
}

