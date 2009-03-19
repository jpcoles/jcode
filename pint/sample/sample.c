#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "tipsydefs.h"

int main(int argc,char **argv)
{
    struct dump h;
    struct dark_particle *dp,dt;
    struct gas_particle *gp,gt;
    struct star_particle *sp,st;
    int i,k,iDum,nDark,nGas,nStar;
    double di;

    iDum = -1;
    i = 1;
    nDark = 0;
    nStar = 0;
    nGas = 0;
    while (i < argc) {
        if (!strcmp(argv[i],"-d")) {
            ++i;
            nDark = atoi(argv[i]);
            ++i;
            }
        else if (!strcmp(argv[i],"-g")) {
            ++i;
            nGas = atoi(argv[i]);
            ++i;
            }
        else if (!strcmp(argv[i],"-s")) {
            ++i;
            nStar = atoi(argv[i]);
            ++i;
            }
        else ++i;
        }
    fread(&h,sizeof(struct dump),1,stdin);
    if (nDark > h.ndark) nDark = h.ndark;
    if (nGas > h.nsph) nGas = h.nsph;
    if (nStar > h.nstar) nStar = h.nstar;

    fprintf(stderr, "nDark=%i nGas=%i nStart=%i\n", nDark, nGas, nStar);

    if (nGas > 0) {
        gp = (struct gas_particle *)malloc(h.nsph*
                           sizeof(struct gas_particle));
        assert(gp != NULL);
        fread(gp,sizeof(struct gas_particle),h.nsph,stdin);
        for (i=h.nsph-1;i > 0;--i) {
            di = i;
            k=(int)(di*rand()/(RAND_MAX+1.0));
            gt = gp[k];
            gp[k] = gp[i];
            gp[i] = gt;
            }
        }
    if (nDark > 0) {
        dp = (struct dark_particle *)malloc(h.ndark*
                            sizeof(struct dark_particle));
        assert(dp != NULL);
        fread(dp,sizeof(struct dark_particle),h.ndark,stdin);
        for (i=h.ndark-1;i > 0;--i) {
            di = i;
            k=(int)(di*rand()/(RAND_MAX+1.0));
            dt = dp[k];
            dp[k] = dp[i];
            dp[i] = dt;
            }
        }
    if (nStar > 0) {
        sp = (struct star_particle *)malloc(h.nstar*
                            sizeof(struct star_particle));
        assert(sp != NULL);
        fread(sp,sizeof(struct star_particle),h.nstar,stdin);
        for (i=h.nstar-1;i > 0;--i) {
            di = i;
            k=(int)(di*rand()/(RAND_MAX+1.0));
            st = sp[k];
            sp[k] = sp[i];
            sp[i] = st;
            }
        }
    h.ndark = nDark;
    h.nsph = nGas;
    h.nstar = nStar;
    h.nbodies = nDark+nStar+nGas;
    fwrite(&h,sizeof(struct dump),1,stdout);
    if (nGas > 0) {
        fwrite(gp,sizeof(struct gas_particle),nGas,stdout);
        free(gp);
        }
    if (nDark > 0) {
        fwrite(dp,sizeof(struct dark_particle),nDark,stdout);
        free(dp);
        }
    if (nStar > 0) {
        fwrite(sp,sizeof(struct star_particle),nStar,stdout);
        free(sp);
        }
    }

