#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include "TP1.h"
#include "SolveurMC.h"

using namespace std;

int main()
{

    //Initialistaion des constantes

    REAL st=1;
    int N=100000;
    int K=100;
    REAL mu=0.5;

    char source[]="Uniforme";

    std::vector<REAL> abs (K+1);
    std::vector<REAL> E (K+1);


    //Appel au solveur

    std::vector<REAL> I=SolveurMC(N, K, st, mu, source);



    //Calcul de la solution exacte

    if (strcmp(source,"Dirac")==0){
      for (int i=0 ; i<abs.size() ; ++i) {
          abs[i] = i/((REAL)K);
          //solution exacte pour source dirac 0
          E[i]=(1/mu)*exp(-st*abs[i]/mu);
        }
    }
    if (strcmp(source,"Uniforme")==0){
      if(mu>=0){
        for (int i=0 ; i<abs.size() ; ++i) {
          abs[i] = i/((REAL)K);
          //solution exacte pour source uniforme, mu positif
          E[i]=(1/st)*(1-exp(-st*abs[i]/mu));
        }
      }
      if(mu<0){
        for (int i=0 ; i<abs.size() ; ++i) {
          abs[i] = i/((REAL)K);
          //solution exacte pour source uniforme, mu négatif
          E[i]=(1/st)*(1-exp(-st*(abs[i]-1)/mu));
        }
      }
    }


     //Ecritures des résultats

    {
        std::ofstream solveur ("solveur.dat");
        std::ofstream exact ("exact.dat");
        for (int i=0 ; i<abs.size() ; ++i) {
        solveur << abs[i] << " " << I[i] << std::endl;
        exact << abs[i] << " " << E[i]<< std::endl;
        }
    }
    //cout<<"u="<<u<<endl;
    //cout<<"u réél="<<2*exp(-2)<<endl;
}
