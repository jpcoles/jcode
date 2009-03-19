#include <math.h>
#include <assert.h>
#include "heidi.h"


extern Environment env;

inline Real SQNORM(Vector(Real, vec))
{
    Real sqmag = 0;
    V( sqmag += vec _ * vec _ );
    return sqmag;
}

#define drift(p, dt) V( (p).x _ += (p).v _ * dt )
#define kick(p, dt)  V( (p).v _ += (p).a _ * dt )
//#define drift(p, dt) V( (p).x _ += (p).v _ * (p).dt )
//#define kick(p, dt)  V( (p).v _ += (p).a _ * (p).dt )

//Real close_approach(Vector(Real, 

void update_accelerations()
{
    Integer i, j;
    Real f, r;

    Vector(Real, dr);

    for all_particles(i) 
        V( P_a(i)_ = 0 );

    for all_particle_pairs(i, j)
    {
        V( dr _ = P_x(i)_ - P_x(j)_ );

        Dist d0 = sqrt(SQNORM(dr)); assert(d0 != 0);
        Dist d1;
        
        r =  1 / d0; //assert(isfinite(r));

        /* Gravity */
        f = _G * P_m(j);
        V( P_a(i)_ += (((f * dr _) * r) * r) * r );

        f = _G * P_m(i);
        V( P_a(j)_ -= (((f * dr _) * r) * r) * r );

        /* Charge */
        f = _kC * _fsc * P_q(i) * P_q(j);
        //fprintf(env.io.err, "f= %e %e %e\n", f, P_q(i), P_q(j));
        V( P_a(i) _ -= (((f * dr _) * r) * r) * r );
        V( P_a(j) _ += (((f * dr _) * r) * r) * r );

        // http://mathworld.wolfram.com/Line-LineDistance.html

        Time dt = MIN(P_dt(i), P_dt(j));

        Real D;
        Vector(Real, a);
        Vector(Real, b);
        Vector(Real, c);
        Vector(Real, X);

        V( a _ = P_v(i) _ * dt );
        V( b _ = P_v(j) _ * dt );
        V( c _ = P_x(j) _ - P_x(i) _ );


        //fprintf(env.io.out, "A = %e %e %e %e\n", A[0], A[1], A[2], NORM(A));
        //fprintf(env.io.out, "B = %e %e %e %e\n", B[0], B[1], B[2], NORM(B));
        //fprintf(env.io.out, "D = %e\n", D);
        //fprintf(env.io.out, "C = %e %e %e %e\n", C[0], C[1], C[2], NORM(C));
        //fprintf(env.io.out, "C = %e %e %e %e\n", Chat[0], Chat[1], Chat[2], NORM(Chat));
        //fprintf(env.io.out, "dot = %e\n", DOT(a, c0));


        VECIF(a _ == 0)
            VECIF(b _ == 0)
                continue;
            ENDVECIF

            V( SWAP(a _, b _) );

            Vector(Real, chat);
            D = DOT(a, c) / NORM(c);
            UNIT(c, chat);
            V( chat _ = chat _ * D );
            V( chat _ = a _ - chat _ );
            UNIT(chat, b);
        ELSEVECIF 
            VECIF(b _ == 0)

                Vector(Real, chat);
                D = DOT(a, c) / NORM(c);
                UNIT(c, chat);
                V( chat _ = chat _ * D );
                V( chat _ = a _ - chat _ );
                UNIT(chat, b);
            ENDVECIF
        ENDVECIF

        CROSS(a, b, X);

        d1 = d0;
        if (NORM(X) != 0)
        {
            d1  = ABS(DOT(c, X));
            d1 /= NORM(X);
        }
        
        //Real r = d0 * d0 / d1 / d1;
        Real r = pow(d0 / d1, 0.25) * SECONDS(2e-13);

        //fprintf(env.io.out, "++++++++++++++++\n");
        //fprintf(env.io.out, "P_x = %e %e %e [dt = %e]\n", P_x(i)[0], P_x(i)[1], P_x(i)[2], P_dt(i));
        //fprintf(env.io.out, "a = %e %e %e %e\n", a[0], a[1], a[2], NORM(a));
        //fprintf(env.io.out, "b = %e %e %e %e\n", b[0], b[1], b[2], NORM(b));
        //fprintf(env.io.out, "c = %e %e %e %e\n", c[0], c[1], c[2], NORM(c));
        //fprintf(env.io.out, "X = %e %e %e %e\n", X[0], X[1], X[2], NORM(X));
        //fprintf(env.io.out, "d = %e\n", d);

        //P_dt(i) = MAX(MIN(P_dt(i), r), SECONDS(2e-50));
        //P_dt(j) = MAX(MIN(P_dt(j), r), SECONDS(2e-50));
        P_dt(i) = r; //MIN(P_dt(i), r);
        P_dt(j) = r; //MIN(P_dt(j), r);

        fprintf(env.io.out, "%e %e %e %e %e %e %e\n", P_xx(0), d0, d1, r, P_dt(i), P_dt(j), SECONDS(2e-10));
    }

    for all_particles(i) 
        V( P_a(i)_ /= P_m(i) );
}

void update_energy()
{
    int i;

    env.E.kin = 0;

    for all_particles(i) 
    {
        Vel v2 = SQNORM(P_v(i)); assert(v2 != 0);

        env.E.kin += 0.5 * P_m(i) * v2;
    }

    env.E.kin /= NP;
}

Time step_particles()
{
    int i;
    Time dt = SECONDS(2e-10);

    for all_particles(i) drift(P(i), dt);
    update_accelerations();
    for all_particles(i) kick(P(i),  dt);
    for all_particles(i) drift(P(i), dt);

    update_energy();

    return dt;
}

