//
// File Manager
//
// Joachim Stadel, 2003.
//
//
// This program and source code are under the
// GNU Public License (GPL)

#ifndef TIPSY_INCLUDED
#define TIPSY_INCLUDED

#include <rpc/types.h>
#include <rpc/xdr.h>

#define TIPSY_TYPE_GAS	1
#define TIPSY_TYPE_DARK	2
#define TIPSY_TYPE_STAR	3

#define TIPSY_BLOCK 10000

struct gas_particle {
    float mass;
    float pos[3];
    float vel[3];
    float rho;
    float temp;
    float hsmooth;
    float metals;
    float phi;
};

struct dark_particle {
    float mass;
    float pos[3];
    float vel[3];
    float eps;
    float phi;
};

struct star_particle {
    float mass;
    float pos[3];
    float vel[3];
    float metals;
    float tform;
    float eps;
    float phi;
};

struct base_particle {
    float mass;
    float pos[3];
    float vel[3];
};

struct dump {
    double time;
    int nbodies;
    int ndim;
    int nsph;
    int ndark;
    int nstar;
};

struct TipsyContext {
	int bNative;
	int iCurr;
	int iCurrType;
	int nRemain;
	int nBlock;
	int iBlock;
	struct gas_particle gpCurr;
	struct dark_particle dpCurr;
	struct star_particle spCurr;
	int nMaxGas,nGas;
	int nMaxDark,nDark;
	int nMaxStar,nStar;
	struct gas_particle *gp;
	struct dark_particle *dp;
	struct star_particle *sp;
	FILE *fp;
	XDR xdr;
	struct dump head;
	char *pchInFile;
};

typedef struct TipsyContext* TCTX;

int xdr_header(XDR *xdrs, struct dump *header);
int xdr_gas(XDR *xdrs,struct gas_particle *p);
int xdr_dark(XDR *xdrs,struct dark_particle *p);
int xdr_star(XDR *xdrs,struct star_particle *p);

int TipsyInitialize(TCTX *pctx,int bNative,char *pchInFile);
void TipsyFinish(TCTX ctx);
struct base_particle *pTipsyReadNative(TCTX ctx,int *piType,double *pdSoft);
struct base_particle *pTipsyRead(TCTX ctx,int *piType,double *pdSoft);
struct base_particle *pTipsyParticle(TCTX ctx,int iIndex,int *piType,double *pdSoft);
int TipsyReadAll(TCTX ctx);
void TipsyAddGas(TCTX ctx,struct gas_particle *pGas);
void TipsyAddDark(TCTX ctx,struct dark_particle *pDark);
void TipsyAddStar(TCTX ctx,struct star_particle *pStar);
void TipsyWriteAll(TCTX ctx,double dTime,char *pchFileName);
int iTipsyNumParticles(TCTX ctx);
double dTipsyTime(TCTX ctx);

#endif

