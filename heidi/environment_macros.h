#ifndef ENVIRONMENT_MACROS_H
#define ENVIRONMENT_MACROS_H

#define _G       (env.K.G)
#define _kC      (env.K.kC)
#define _kB      (env.K.kB)
#define _fsc     (env.K.fsc)

#define _Me      (env.K.Me)
#define _Mp      (env.K.Mp)
#define _Mn      (env.K.Mn)

#define _Qe      (env.K.Qe)
#define _Qp      (env.K.Qp)
#define _Qn      (env.K.Qn)

#define _h       (env.K.h)
#define _hbar    (env.K.hbar)

#define _c       (env.K.c)


#define NP       (env.plist.N)
#define P(i)     (env.plist.p[i])

#define P_type(i) (P(i).type)
#define P_x(i)   (P(i).x)
#define P_v(i)   (P(i).v)
#define P_a(i)   (P(i).a)
#define P_dt(i)  (P(i).dt)
#define P_m(i)   (env.particle_type_data[P_type(i)].m)
#define P_q(i)   (env.particle_type_data[P_type(i)].q)
#define P_Z(i)   (env.particle_type_data[P_type(i)].Z)
#define P_N(i)   (env.particle_type_data[P_type(i)].N)
#define P_E(i)   (env.particle_type_data[P_type(i)].E)

#define P_xx(i)  (VECELEM(P_x(i), 0))
#define P_xy(i)  (VECELEM(P_x(i), 1))
#define P_xz(i)  (VECELEM(P_x(i), 2))
#define P_vx(i)  (VECELEM(P_v(i), 0))
#define P_vy(i)  (VECELEM(P_v(i), 1))
#define P_vz(i)  (VECELEM(P_v(i), 2))
#define P_ax(i)  (VECELEM(P_a(i), 0))
#define P_ay(i)  (VECELEM(P_a(i), 1))
#define P_az(i)  (VECELEM(P_a(i), 2))

#endif
