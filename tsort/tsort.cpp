#include <iostream.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <ftipsy.hpp>

void help()
{
    cerr << "Usage: tsort <tipsy-standard-file>" << endl;
    exit(2);
}

int compar(const void *a0, const void *b0)
{
    TipsyBaseParticle *a = (TipsyBaseParticle*)a0;
    TipsyBaseParticle *b = (TipsyBaseParticle*)b0;

    float ar = a->pos[0]*a->pos[0] + a->pos[1]*a->pos[1] + a->pos[2]*a->pos[2];
    float br = b->pos[0]*b->pos[0] + b->pos[1]*b->pos[1] + b->pos[2]*b->pos[2];

    if (ar > br) return 1;
    if (ar < br) return 1;

    return 0;
}

int main(int argc, char **argv)
{
    if (argc < 2) help();

    ifTipsy in;

    TipsyHeader h;
    TipsyGasParticle g;
    TipsyDarkParticle d;
    TipsyStarParticle s;

    in.open(argv[1], "standard");

    if (!in.is_open())
    {
        cerr << "Can't open file " << argv[1] << endl;
        exit(1);
    }

    in >> h;
    
    TipsyBaseParticle **ps = new TipsyBaseParticle*[h.h_nBodies];
    assert(ps != NULL);

    unsigned int j=0;
    for (unsigned int i=0; i < h.h_nGas;  i++) { ps[j] = new TipsyGasParticle();  in >> g; *(ps[j]) = g; j++; }
    for (unsigned int i=0; i < h.h_nDark; i++) { ps[j] = new TipsyDarkParticle(); in >> d; *(ps[j]) = d; j++; }
    for (unsigned int i=0; i < h.h_nStar; i++) { ps[j] = new TipsyStarParticle(); in >> s; *(ps[j]) = s; j++; }

    in.close();

    //qsort(ps, h.h_nBodies, sizeof(TipsyBaseParticle*), compar);

    unsigned int k=0;
    float EPS = 1e-1;
    for (unsigned int i=0; i < h.h_nBodies-1; i++)
    {
        for (unsigned int j=i+1; j < h.h_nBodies; j++)
        {
            if (fabsf(ps[i]->pos[0] - ps[j]->pos[0]) < EPS && 
                fabsf(ps[i]->pos[1] - ps[j]->pos[1]) < EPS && 
                fabsf(ps[i]->pos[2] - ps[j]->pos[2]) < EPS)
            {
                cout << ps[i]->pos[0] << " " << ps[i]->pos[1] << " " << ps[i]->pos[2] << endl;
                cout << ps[j]->pos[0] << " " << ps[j]->pos[1] << " " << ps[j]->pos[2] << endl;
                k++;
            }
        }
    }

    cout << k << " particles are close to one another." << endl;

    return 0;
}
