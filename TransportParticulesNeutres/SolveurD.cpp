#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include "SolveurD.h"


using namespace std;

std::vector<REAL> SchemaDiamant(int K, REAL mu,std::vector<REAL> sd, std::vector<REAL> Qt,REAL CL){

  std::vector<REAL> PHI(K+1);

  REAL DELTA=1/((REAL)K);

  if (mu>0){
    PHI[0]=CL;
    for (int i=1; i<=K; i++){
      PHI[i]= (2*DELTA*Qt[i-1] + (2*mu - DELTA*sd[i])*PHI[i-1])/(2*mu + DELTA*sd[i]);
    }
  }

  if (mu<0){
    PHI[K]=CL;
    for (int i=K-1; i>=0; i--){
      PHI[i]=(2*DELTA*Qt[i] - (2*mu + DELTA*sd[i])*PHI[i+1])/(DELTA*sd[i] - 2*mu);
    }
  }

  if (mu==0){
    std::cerr << "erreur : mu =0" <<endl;//erreur
  }
 return PHI;

}



std::vector<REAL> MethodeSN(int K,REAL epsilon, REAL CL, std::vector<REAL> sd, std::vector<REAL> S){

  int i=0;
  int j=0;
  int k=0;
  int k2=0;
  std::vector<REAL> phi(K+1);
  std::vector<REAL> phim(K);
  std::vector<REAL> Qt(K);
  std::vector<REAL> Qt_aux(K);
  std::vector<REAL> Qtp(K);
  REAL Qtdiff;

  //calcul des poids de la formule de qudrature de Gauss-Legendre
  REAL Nmu=25;//taille de la discrétisation pour mu, Paire
  std::vector<REAL> mu(2*Nmu);//discrétisation de mu
  for (int i=0; i<2*Nmu; i++){
    mu[i]=-1+1/(2*Nmu)+i/Nmu;
    cout<<"mu["<<i<<"]="<<mu[i]<<endl;
  }

  REAL P=1/(2*Nmu);

  Qt=S;
  Qt_aux=S;


  while ((k==0)||(k2==0)||(Qtdiff>epsilon)){
    Qtdiff=0;

    if (k==1){
        k2=1;//pour rentrer dans la boucle a la deuxieme itération
    }
    k=1;//pour rentrer dans le while à la premere itération

    for (i=0; i<2*Nmu; i++){
       phi=SchemaDiamant( K, mu[i], sd, Qt, CL);
       for (int j=0; j<K; j++){
         phim[j]=(phi[j]+phi[j+1])/2;
         Qt_aux[j]=Qt_aux[j]+sd[j]*P*phim[j];
       }
       //equation dsa Qt_aux=Qt_aux + ......
    }
    Qtp=Qt;
    Qt=Qt_aux;
    Qt_aux=S;

    for (j=0; j<K; j++){
      Qtdiff=Qtdiff + pow(Qt[j]-Qtp[j],2);
    }
    Qtdiff=sqrt(Qtdiff);

    cout<<"k="<<k<<";k2="<<k2<<";Qtdiff="<<Qtdiff<<";epsilon="<<epsilon<<endl;
  }

  for (j=0; j<K; j++){
    phim[j]=(Qt[j]-S[j])/sd[j];
  }
  return phim;
}
