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
      PHI[i]=(2*DELTA*Qt[i+1] - (2*mu + DELTA*sd[i])*PHI[i+1])/(DELTA*sd[i] - 2*mu);
    }
  }

  if (mu==0){
    std::cerr << "erreur : mu =0" <<endl;//erreur
  }
 return PHI;

}

