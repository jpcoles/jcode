#include <iostream>
#include <vector>
#include <stdlib.h>
#include "ftipsy.hpp"

using namespace std;

#define N 10

class Density
{
public:
    double den;
    uint64_t index;

    Density(double d, uint64_t i) : den(d), index(i) {}
};

void help()
{
    cerr << "Usage: denminmax2 <tipsy> <density>" << endl;
    exit(2);
}

int main(int argc, char **argv)
{
    if (argc < 3) help();

    string tipsyfile(argv[1]);
    string denfile(argv[2]);

    uint64_t nParticles;
    double den;
    //double minden = 1e30, maxden = 0;
    //uint64_t min_index=-1, max_index=-1;

    vector<Density> minden;
    vector<Density> maxden;

    ifstream in;
    in.open(denfile.c_str(), ios::in);
    in >> nParticles;
    for (uint64_t i=0; in >> den; i++)
    {
        if (1e-7 <= den&&den <= 1e-6) continue;

        int best=-1;
        int maxi=0;
        for (uint32_t j=0; j < minden.size(); j++)
        {
            if (den  < minden[j].den) best = j;
            if (maxi > minden[j].den) maxi = j;
        }

        if (minden.size() < N) minden.insert(minden.begin() + (best==-1?0:best), Density(den, i));
        else if (best != -1) minden[best] = Density(den, i);

        best=-1;
        int mini=0;
        for (uint32_t j=0; j < maxden.size(); j++)
        {
            if (den  > maxden[j].den) best = j;
            if (mini < maxden[j].den) mini = j;
        }

        if (maxden.size() < N) maxden.insert(maxden.begin() + (best==-1?0:best), Density(den, i));
        else if (best != -1) maxden[best] = Density(den, i);

        //if (den < minden) { minden = den; min_index = i; }
        //if (den > maxden) { maxden = den; max_index = i; }
    }
    in.close();

    if (minden.size() == 0 || maxden.size() == 0)
    {
        cerr << "No min/max found. Corrupt density file?" << endl;
        exit(1);
    }

    ifTipsy tin;
    TipsyHeader h;
    TipsyDarkParticle d;

    tin.open(tipsyfile.c_str());

    tin >> h;

    cout << "Minimum density " << endl;
    for (uint32_t j=0; j < minden.size(); j++)
    {
        tin.seekg(tipsypos(tipsypos::particle, minden[j].index));
        tin >> d;
        cout << "\t" 
             << minden[j].den << " at particle " << minden[j].index 
             << " ( " << d.pos[0] << " " << d.pos[1] << " " << d.pos[2] << " )" << endl;
    }

    double r=0.01;
    for (uint32_t j=0; j < minden.size(); j++)
    {
        tin.seekg(tipsypos(tipsypos::particle, minden[j].index));
        tin >> d;
        cout << "\tsetbox " << (j+1) << " " << d.pos[0] << " " << d.pos[1] << " " << d.pos[2] << " ";
        cout << r << " " << r << " " << r << endl;
        r += 0.01;
    }

    cout << endl;
    cout << "Maximum density " << endl;
    for (uint32_t j=0; j < maxden.size(); j++)
    {
        tin.seekg(tipsypos(tipsypos::particle, maxden[j].index));
        tin >> d;
        cout << "\t" 
             << maxden[j].den << " at particle " << maxden[j].index 
             << " ( " << d.pos[0] << " " << d.pos[1] << " " << d.pos[2] << " )" << endl;
    }

    r=0.01;
    for (uint32_t j=0; j < maxden.size(); j++)
    {
        tin.seekg(tipsypos(tipsypos::particle, maxden[j].index));
        tin >> d;
        cout << "\tsetbox " << (j+1) << " " << d.pos[0] << " " << d.pos[1] << " " << d.pos[2] << " ";
        cout << r << " " << r << " " << r << endl;
        r += 0.01;
    }

    tin.close();

    return 0;
}

