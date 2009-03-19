#include <stdlib.h>
#include <assert.h>
#include <string>
#include <math.h>
#include <complex>
#include <iostream>
#include <pthread.h>
#include "ftipsy.hpp"

using namespace std;

#if 0

#if 1
#define NTHREADS 10

#define K ((int)8)
#define L ((double)1.0)

typedef struct
{
    uint32_t id;
    complex<double> C[2*K+1][2*K+1][2*K+1];
    string file, mode;
    uint64_t offs, len;
} worker_t;


void help()
{
    cerr << "Usage: denminmax <tipsy-in>" << endl;
    exit(2);
}

void *worker(void *arg)
{
    worker_t *w = (worker_t *)arg;

    const complex<double> IPI_L(0, M_PI/L);

    ifTipsy in;
    TipsyHeader h;
    TipsyDarkParticle d;

    in.open(w->file.c_str(), w->mode.c_str());
    if (!in.is_open())
    {
        cerr << "Can't open input file " << w->file << endl;
        exit(1);
    }

    in >> h;
    in.seekg(tipsypos(tipsypos::particle, w->offs));
    for (int i=0; i < w->len; i++)
    {
        in >> d;
        float x = d.pos[0];
        float y = d.pos[1];
        float z = d.pos[2];

        for (int kz=-K; kz <= K; kz++)
        {
            double kzz = kz*z;
            for (int ky=-K; ky <= K; ky++)
            {
                double kyy = ky*y;
                for (int kx=-K; kx <= K; kx++)
                {
                    double k_dot_x = kx*x + kyy + kzz;
                    w->C[kz+K][ky+K][kx+K] += exp(-IPI_L * k_dot_x);
                }
            }
        }

        if ((i%100000) == 0) cerr << w->id;
        //if ((i%10000) == 0) cerr << ((100.0 * i) / h.h_nDark) << endl;
    }

    cerr << endl << "Thread " << w->id << " finished." << endl;

    in.close();

    pthread_exit(NULL);
}
#endif

int main(int argc, char **argv)
{

#if 1

    ifTipsy in;
    TipsyHeader h;
    TipsyDarkParticle d;
    TipsyStarParticle s;
    TipsyGasParticle  g;
    string mode = "standard";

    if (argc < 2) help();

    string input_filename(argv[1]);

    in.open(input_filename.c_str(), mode.c_str());
    if (!in.is_open())
    {
        cerr << "Can't open input file " << input_filename << endl;
        exit(1);
    }

    in >> h;
    in.close();

    pthread_t *thrds = new pthread_t[NTHREADS];
    worker_t  *ws    = new worker_t[NTHREADS];
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    uint64_t offs = 0;
    uint64_t len  = min(h.h_nDark / NTHREADS + 1, h.h_nDark);

    cerr << "nBodies=" << h.h_nBodies << endl;

    for (int i=0; i < NTHREADS; i++)
    {
        ws[i].id   = i;
        ws[i].file = input_filename;
        ws[i].mode = mode;
        ws[i].offs = offs;
        ws[i].len  = (offs + len < h.h_nDark) ? len : h.h_nDark - offs;
        cerr << i << "] offs=" << ws[i].offs 
                  << " len=" << ws[i].len 
                  << " sum=" << (ws[i].offs+ws[i].len) << endl;
        memset(ws[i].C, 0, sizeof(ws[i].C));
        offs += len;
        if (pthread_create(&thrds[i], &attr, worker, (void *)&ws[i]) != 0)
        {
            cerr << "Error creating threads." << endl;
            exit(1);
        }
    }

    complex<double> C[2*K+1][2*K+1][2*K+1];
    memset(C, 0, sizeof(C));

    void *status;
    pthread_attr_destroy(&attr);
    for (int i=0; i < NTHREADS; i++)
    {
        pthread_join(thrds[i], &status);
        for (int kz=-K; kz <= K; kz++)
            for (int ky=-K; ky <= K; ky++)
                for (int kx=-K; kx <= K; kx++)
                    C[kz+K][ky+K][kx+K] += ws[i].C[kz+K][ky+K][kx+K];

    }

#else

#   include "C6"

#endif
    cout << "const int     K = " << K << ";" << endl;
    cout << "const double  L = " << L << ";" << endl;
    cout << "const complex<double> IPI_L(0, M_PI/L);" << endl;
    cout << "const complex<double> C[2*K+1][2*K+1][2*K+1] = " << endl;
    cout << "{ ";
    for (int kz=-K; kz <= K; kz++)
    {
        cout << "{ ";
        for (int ky=-K; ky <= K; ky++)
        {
            cout << "{ ";
            for (int kx=-K; kx <= K; kx++)
            {
                //cout << "C[" << (kz+K) << "][" << (ky+K) << "][" << (kx+K) << "] = ";
                cout << "complex<double>("
                     << C[kz+K][ky+K][kx+K].real() << ", " 
                     << C[kz+K][ky+K][kx+K].imag() << ")";
                cout << "," << endl;
                //cout << ";" << endl;
            }
            cout << "}, ";
        }
        cout << "}, ";
    }
    cout << "};";

    return 0;
}

#else

int main(int argc, char **argv)
{
#   include "C7"
#   define GRID .01

    double min[4] = {0,0,0,1e30}, max[4] = {0,0,0,-1e30};

    for (double z=-L/2; z <= L/2; z += GRID)
    {
        for (double y=-L/2; y <= L/2; y += GRID)
        {
            for (double x=-L/2; x <= L/2; x += GRID)
            {
                double p=0;

                for (int kz=-K; kz <= K; kz++)
                {
                    double kzz = kz*z;
                    for (int ky=-K; ky <= K; ky++)
                    {
                        double kyy = ky*y;
                        for (int kx=-K; kx <= K; kx++)
                        {
                            double k_dot_x = kx*x + kyy + kzz;
                            p += (C[kz+K][ky+K][kx+K] * exp(IPI_L * k_dot_x)).real();
                        }
                    }
                }

                if (p < min[3]) 
                {
                    min[0] = x;
                    min[1] = y;
                    min[2] = z;
                    min[3] = p;
                }

                if (p > max[3]) 
                {
                    max[0] = x;
                    max[1] = y;
                    max[2] = z;
                    max[3] = p;
                }
            }
            cerr << ".";
        }
        cerr << endl << "z=" << z << "|";
    }
    cerr << endl;

    printf("Min density %e at (%f %f %f).\n", min[3], min[0],min[1],min[2]);
    printf("Max density %e at (%f %f %f).\n", max[3], max[0],max[1],max[2]);
    return 0;
}

#endif

