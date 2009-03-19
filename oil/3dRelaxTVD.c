#include <stdio.h>
#include <math.h>
#include <string.h>
#include "timing.h"

#include "3dRelaxTVD.h"

CPUDEFS

FILE *in, *out, *err;

real gamma = 5.0F/3;
real cfl = 0.9F;

int nsw;
int stopsim=0;
real t;
real tf;
real dt;
real E0;
real rmax;

cell_t u[nc][nc][nc];
cell_t u1d[nc];

limiter_fn_t limiter;

static int ic_sedovtaylor()
{
    int x,y,z;

    for (z=0; z < nc; z++)
    {
        for (y=0; y < nc; y++)
        {
            for (x=0; x < nc; x++)
            {
                u[z][y][x].p   = 1;
                u[z][y][x].pvx = 0;
                u[z][y][x].pvy = 0;
                u[z][y][x].pvz = 0;
                u[z][y][x].e   = 1e-3F;
            }
        }
    }

    u[hc-1][hc-1][hc-1].e = E0;

    return 0;
}

static void write_output()
{
}

static int timestep()
{
    int x,y,z;

    real P, cs, cmax;
    real v[3];

    cmax = 1e-5F;

    for (z=0; z < nc; z++)
    {
        for (y=0; y < nc; y++)
        {
            for (x=0; x < nc; x++)
            {
                v[0] = abs(u[z][y][x].pvx / u[z][y][x].p);
                v[1] = abs(u[z][y][x].pvy / u[z][y][x].p);
                v[2] = abs(u[z][y][x].pvz / u[z][y][x].p);

                real sum_of_sqs = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
                P = max((gamma-1) * (u[z][y][x].e - u[z][y][x].p * sum_of_sqs / 2), 0.F);
                cs = sqrt(gamma * P / u[z][y][x].p);
                cmax = max(cmax, max(max(v[0], v[1]), v[2]) + cs);
            }
        }
    }

    dt = cfl / cmax;
    if (t + 2*dt > tf) 
    {
        dt = (tf-t) / 2;
        stopsim = 1;
    }

    t += 2 * dt;
    nsw++;

    fprintf(out, "nsw=%i dt=%f t=%f\n", nsw, dt, t);

    return 0;
}

static int averageflux(cell_t *u, cell_t *w, real *c)
{
    int i;
    real P, v, sum;

    for (i=0; i < nc; i++)
    {
        v = u[i].pvx / u[i].p;
        sum = u[i].pvx * u[i].pvx 
            + u[i].pvy * u[i].pvy 
            + u[i].pvz * u[i].pvz;
        P = max((gamma-1) * (u[i].e - sum / u[i].p / 2), 0.0F);

        c[i] = abs(v) + max(sqrt(gamma * P / u[i].p), 1e-5F);
        w[i].p   = (u[i].p   * v);
        w[i].pvx = (u[i].pvx * v) + P;
        w[i].pvy = (u[i].pvy * v);
        w[i].pvz = (u[i].pvz * v);
        w[i].e   = (u[i].e   + P) * v;
    }

    return 0;
}

static int relaxing(cell_t *u)
{
    int i;
    real c[nc];
    cell_t u1[nc];
    cell_t w[nc];
    cell_t fu[nc];
    cell_t fr[nc];
    cell_t fl[nc];
    cell_t dfl[nc];
    cell_t dfr[nc];

    /***********************************************************************/

    averageflux(u, w, c);

    /***********************************************************************/
    for (i=0; i < nc-1; i++)
    {
        fu[i].p   = ((u[i].p   * c[i] + w[i].p)   - (u[i+1].p   * c[i+1] - w[i+1].p)  ) / 2;
        fu[i].pvx = ((u[i].pvx * c[i] + w[i].pvx) - (u[i+1].pvx * c[i+1] - w[i+1].pvx)) / 2;
        fu[i].pvy = ((u[i].pvy * c[i] + w[i].pvy) - (u[i+1].pvy * c[i+1] - w[i+1].pvy)) / 2;
        fu[i].pvz = ((u[i].pvz * c[i] + w[i].pvz) - (u[i+1].pvz * c[i+1] - w[i+1].pvz)) / 2;
        fu[i].e   = ((u[i].e   * c[i] + w[i].e)   - (u[i+1].e   * c[i+1] - w[i+1].e)  ) / 2;
    }
        fu[i].p   = ((u[i].p   * c[i] + w[i].p)   - (u[0].p   * c[0] - w[0].p)  ) / 2;
        fu[i].pvx = ((u[i].pvx * c[i] + w[i].pvx) - (u[0].pvx * c[0] - w[0].pvx)) / 2;
        fu[i].pvy = ((u[i].pvy * c[i] + w[i].pvy) - (u[0].pvy * c[0] - w[0].pvy)) / 2;
        fu[i].pvz = ((u[i].pvz * c[i] + w[i].pvz) - (u[0].pvz * c[0] - w[0].pvz)) / 2;
        fu[i].e   = ((u[i].e   * c[i] + w[i].e)   - (u[0].e   * c[0] - w[0].e)  ) / 2;

    /***********************************************************************/

        u1[0].p   = u[0].p   - (fu[0].p   - fu[nc-1].p)   * dt / 2;
        u1[0].pvx = u[0].pvx - (fu[0].pvx - fu[nc-1].pvx) * dt / 2;
        u1[0].pvy = u[0].pvy - (fu[0].pvy - fu[nc-1].pvy) * dt / 2;
        u1[0].pvz = u[0].pvz - (fu[0].pvz - fu[nc-1].pvz) * dt / 2;
        u1[0].e   = u[0].e   - (fu[0].e   - fu[nc-1].e)   * dt / 2;
    for (i=1; i < nc; i++)
    {
        u1[i].p   = u[i].p   - (fu[i].p   - fu[i-1].p)   * dt / 2;
        u1[i].pvx = u[i].pvx - (fu[i].pvx - fu[i-1].pvx) * dt / 2;
        u1[i].pvy = u[i].pvy - (fu[i].pvy - fu[i-1].pvy) * dt / 2;
        u1[i].pvz = u[i].pvz - (fu[i].pvz - fu[i-1].pvz) * dt / 2;
        u1[i].e   = u[i].e   - (fu[i].e   - fu[i-1].e)   * dt / 2;
    }


    /***********************************************************************/

    averageflux(u1, w, c);

    /***********************************************************************/

    for (i=0; i < nc; i++)
    {
        fr[i].p   = (u1[i].p   * c[i] + w[i].p)   / 2;
        fr[i].pvx = (u1[i].pvx * c[i] + w[i].pvx) / 2;
        fr[i].pvy = (u1[i].pvy * c[i] + w[i].pvy) / 2;
        fr[i].pvz = (u1[i].pvz * c[i] + w[i].pvz) / 2;
        fr[i].e   = (u1[i].e   * c[i] + w[i].e)   / 2;
    }

    /***********************************************************************/

        dfl[0].p   = (fr[0].p   - fr[nc-1].p)   / 2;
        dfl[0].pvx = (fr[0].pvx - fr[nc-1].pvx) / 2;
        dfl[0].pvy = (fr[0].pvy - fr[nc-1].pvy) / 2;
        dfl[0].pvz = (fr[0].pvz - fr[nc-1].pvz) / 2;
        dfl[0].e   = (fr[0].e   - fr[nc-1].e)   / 2;
    for (i=1; i < nc; i++)
    {
        dfl[i].p   = (fr[i].p   - fr[i-1].p)   / 2;
        dfl[i].pvx = (fr[i].pvx - fr[i-1].pvx) / 2;
        dfl[i].pvy = (fr[i].pvy - fr[i-1].pvy) / 2;
        dfl[i].pvz = (fr[i].pvz - fr[i-1].pvz) / 2;
        dfl[i].e   = (fr[i].e   - fr[i-1].e)   / 2;
    }

    /***********************************************************************/

    for (i=0; i < nc-1; i++)
    {
        dfr[i].p   = dfl[i+1].p;
        dfr[i].pvx = dfl[i+1].pvx;
        dfr[i].pvy = dfl[i+1].pvy;
        dfr[i].pvz = dfl[i+1].pvz;
        dfr[i].e   = dfl[i+1].e;
    }
        dfr[nc-1].p   = dfl[0].p;
        dfr[nc-1].pvx = dfl[0].pvx;
        dfr[nc-1].pvy = dfl[0].pvy;
        dfr[nc-1].pvz = dfl[0].pvz;
        dfr[nc-1].e   = dfl[0].e;

    /***********************************************************************/

    limiter(fr, dfl, dfr);

    /***********************************************************************/

    for (i=0; i < nc-1; i++)
    {
        fl[i].p   = (u1[i+1].p   * c[i+1] - w[i+1].p)   / 2;
        fl[i].pvx = (u1[i+1].pvx * c[i+1] - w[i+1].pvx) / 2;
        fl[i].pvy = (u1[i+1].pvy * c[i+1] - w[i+1].pvy) / 2;
        fl[i].pvz = (u1[i+1].pvz * c[i+1] - w[i+1].pvx) / 2;
        fl[i].e   = (u1[i+1].e   * c[i+1] - w[i+1].e)   / 2;
    }
        fl[nc-1].p   = (u1[0].p   * c[0] - w[0].p)   / 2;
        fl[nc-1].pvx = (u1[0].pvx * c[0] - w[0].pvx) / 2;
        fl[nc-1].pvy = (u1[0].pvy * c[0] - w[0].pvy) / 2;
        fl[nc-1].pvz = (u1[0].pvz * c[0] - w[0].pvz) / 2;
        fl[nc-1].e   = (u1[0].e   * c[0] - w[0].e)   / 2;

    /***********************************************************************/

        dfl[0].p   = (fl[nc-1].p   - fl[0].p)   / 2;
        dfl[0].pvx = (fl[nc-1].pvx - fl[0].pvx) / 2;
        dfl[0].pvy = (fl[nc-1].pvy - fl[0].pvy) / 2;
        dfl[0].pvz = (fl[nc-1].pvz - fl[0].pvz) / 2;
        dfl[0].e   = (fl[nc-1].e   - fl[0].e)   / 2;
    for (i=1; i < nc; i++)
    {
        dfl[i].p   = (fl[i-1].p   - fl[i].p)   / 2;
        dfl[i].pvx = (fl[i-1].pvx - fl[i].pvx) / 2;
        dfl[i].pvy = (fl[i-1].pvy - fl[i].pvy) / 2;
        dfl[i].pvz = (fl[i-1].pvz - fl[i].pvz) / 2;
        dfl[i].e   = (fl[i-1].e   - fl[i].e)   / 2;
    }

    /***********************************************************************/

    for (i=0; i < nc-1; i++)
    {
        dfr[i].p   = dfl[i+1].p;
        dfr[i].pvx = dfl[i+1].pvx;
        dfr[i].pvy = dfl[i+1].pvy;
        dfr[i].pvz = dfl[i+1].pvz;
        dfr[i].e   = dfl[i+1].e;
    }
        dfr[nc-1].p   = dfl[0].p;
        dfr[nc-1].pvx = dfl[0].pvx;
        dfr[nc-1].pvy = dfl[0].pvy;
        dfr[nc-1].pvz = dfl[0].pvz;
        dfr[nc-1].e   = dfl[0].e;

    /***********************************************************************/

    limiter(fl, dfl, dfr);

    /***********************************************************************/

    for (i=0; i < nc; i++)
    {
        fu[i].p   = fr[i].p   - fl[i].p;
        fu[i].pvx = fr[i].pvx - fl[i].pvx;
        fu[i].pvy = fr[i].pvy - fl[i].pvy;
        fu[i].pvz = fr[i].pvz - fl[i].pvz;
        fu[i].e   = fr[i].e   - fl[i].e;
    }

    /***********************************************************************/

        u[0].p   -= (fu[0].p   - fu[nc-1].p)   * dt;
        u[0].pvx -= (fu[0].pvx - fu[nc-1].pvx) * dt;
        u[0].pvy -= (fu[0].pvy - fu[nc-1].pvy) * dt;
        u[0].pvz -= (fu[0].pvz - fu[nc-1].pvz) * dt;
        u[0].e   -= (fu[0].e   - fu[nc-1].e)   * dt;
    for (i=1; i < nc; i++)
    {
        u[i].p   -= (fu[i].p   - fu[i-1].p)   * dt;
        u[i].pvx -= (fu[i].pvx - fu[i-1].pvx) * dt;
        u[i].pvy -= (fu[i].pvy - fu[i-1].pvy) * dt;
        u[i].pvz -= (fu[i].pvz - fu[i-1].pvz) * dt;
        u[i].e   -= (fu[i].e   - fu[i-1].e)   * dt;
    }

    /***********************************************************************/

    return 0;
}

static int sweepx()
{
    int y,z;

    for (z=0; z < nc; z++)
    {
        for (y=0; y < nc; y++)
        {
            memcpy(u1d, u[z][y], nc * sizeof(cell_t));

            relaxing(u1d);

            memcpy(u[z][y], u1d, nc * sizeof(cell_t));
        }
    }
    return 0;
}

static int sweepy()
{
    int x,y,z;

    for (z=0; z < nc; z++)
    {
        for (x=0; x < nc; x++)
        {
            for (y=0; y < nc; y++)
            {
                u1d[y].p   = u[z][y][x].p;
                u1d[y].pvx = u[z][y][x].pvy;
                u1d[y].pvy = u[z][y][x].pvx;
                u1d[y].pvz = u[z][y][x].pvz;
                u1d[y].e   = u[z][y][x].e;
            }

            relaxing(u1d);

            for (y=0; y < nc; y++)
            {
                u[z][y][x].p   = u1d[y].p;  
                u[z][y][x].pvx = u1d[y].pvy;
                u[z][y][x].pvy = u1d[y].pvx;
                u[z][y][x].pvz = u1d[y].pvz;
                u[z][y][x].e   = u1d[y].e;
            }
        }
    }

    return 0;
}

static int sweepz()
{
    int x,y,z;

    for (y=0; y < nc; y++)
    {
        for (x=0; x < nc; x++)
        {
            for (z=0; z < nc; z++)
            {
                u1d[z].p   = u[z][y][x].p;
                u1d[z].pvx = u[z][y][x].pvz;
                u1d[z].pvy = u[z][y][x].pvy;
                u1d[z].pvz = u[z][y][x].pvx;
                u1d[z].e   = u[z][y][x].e;
            }

            relaxing(u1d);

            for (z=0; z < nc; z++)
            {
                u[z][y][x].p   = u1d[z].p;  
                u[z][y][x].pvx = u1d[z].pvz;
                u[z][y][x].pvy = u1d[z].pvy;
                u[z][y][x].pvz = u1d[z].pvx;
                u[z][y][x].e   = u1d[z].e;
            }
        }
    }

    return 0;
}


static int vanleer(cell_t *f, cell_t *a, cell_t *b)
{
    int i;
    real c;

    for (i=0; i < nc; i++)
    {
        if ((c = a[i].p   * b[i].p)   > 0) f[i].p   += 2 * c / (a[i].p   + b[i].p);
        if ((c = a[i].pvx * b[i].pvx) > 0) f[i].pvx += 2 * c / (a[i].pvx + b[i].pvx);
        if ((c = a[i].pvy * b[i].pvy) > 0) f[i].pvy += 2 * c / (a[i].pvy + b[i].pvy);
        if ((c = a[i].pvz * b[i].pvz) > 0) f[i].pvz += 2 * c / (a[i].pvz + b[i].pvz);
        if ((c = a[i].e   * b[i].e)   > 0) f[i].e   += 2 * c / (a[i].e   + b[i].e);
    }

    return 0;
}

#define SIGN(a) ((a) > 0 ? 1 : 0)
#define MINMOD(a,b, el) \
    ((SIGN(a[i]. el) + SIGN(b[i]. el)) * min(abs(a[i]. el), abs(b[i]. el)) / 2)

static int minmod(cell_t *f, cell_t *a, cell_t *b)
{
    int i;

    for (i=0; i < nc; i++)
    {
        f[i].p   += MINMOD(a, b, p);
        f[i].pvx += MINMOD(a, b, pvx);
        f[i].pvy += MINMOD(a, b, pvy);
        f[i].pvz += MINMOD(a, b, pvz);
        f[i].e   += MINMOD(a, b, e);
    }

    return 0;
}

static int superbee(cell_t *f, cell_t *a, cell_t *b)
{
    int i;

    for (i=0; i < nc; i++)
    {
        f[i].p   += (abs(a[i].p)   >= abs(b[i].p))   ? MINMOD(a, 2*b, p)   : MINMOD(2*a, b, p);
        f[i].pvx += (abs(a[i].pvx) >= abs(b[i].pvx)) ? MINMOD(a, 2*b, pvx) : MINMOD(2*a, b, pvx);
        f[i].pvy += (abs(a[i].pvy) >= abs(b[i].pvy)) ? MINMOD(a, 2*b, pvy) : MINMOD(2*a, b, pvy);
        f[i].pvz += (abs(a[i].pvz) >= abs(b[i].pvz)) ? MINMOD(a, 2*b, pvz) : MINMOD(2*a, b, pvz);
        f[i].e   += (abs(a[i].e)   >= abs(b[i].e))   ? MINMOD(a, 2*b, e)   : MINMOD(2*a, b, e);
    }

    return 0;
}

int main(int argc, char **argv)
{
    double start_time, end_time;

    in  = stdin;
    out = stdout;
    err = stderr;

    t = 0;
    dt = 0;
    nsw = 0;
    stopsim = 0;

    E0 = 1e5F;
    rmax = 3*hc / 4;
    tf = sqrt(pow(rmax / 1.15F, 5) / E0);
    
    fprintf(out, "tf = %f\n", tf);

    limiter = vanleer;
    
    ic_sedovtaylor();

    start_time = CPUTIME;

    while (1)
    {
        timestep();
        sweepx();
        sweepy();
        sweepz();
        sweepz();
        sweepy();
        sweepx();

        if (stopsim) break;

        timestep();
        sweepz();
        sweepx();
        sweepy();
        sweepy();
        sweepx();
        sweepz();

        if (stopsim) break;

        timestep();
        sweepy();
        sweepz();
        sweepx();
        sweepx();
        sweepz();
        sweepy();

        if (stopsim) break;
    }

    end_time = CPUTIME;

    FILE *time_fp = fopen("3dRelaxTVD.timing", "a");
    fprintf(time_fp, "%f\n", (end_time - start_time));
    fclose(time_fp);

    write_output();

    return 0;
}
