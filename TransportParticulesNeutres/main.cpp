#include <iostream>
#include <fstream>
#include <math.h>
#include "TP1.h"
#include "SolveurMC.h"

using namespace std;

int main()
{
    REAL st=1;
    int N=10000;
    int K=100;

    std::vector<REAL> I=SolveurMC(N, K, st);

    std::vector<REAL> abs (K+1);
    for (int i=0 ; i<abs.size() ; ++i) {
        abs[i] = i/((REAL)K);
    }


    {
        std::ofstream toto ("monfichier.dat");
        for (int i=0 ; i<abs.size() ; ++i) {
        toto << abs[i] << " " << I[i] << std::endl;
        }
    }
    //cout<<"u="<<u<<endl;
    //cout<<"u réél="<<2*exp(-2)<<endl;
}
