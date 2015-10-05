#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include<vector>
#include "SolveurMC.h"


using namespace std;

std::vector<REAL> SolveurMC(int N, int K, REAL st){

  std::vector<REAL> X (N);
  std::vector<REAL> I(K+1);
  REAL mu=0.5;
  REAL Q=0;
  REAL u=0;
  int j;

  srand (time(NULL));

//initialisation de I

  for (int i=0; i<K; i=i+1){
    I[i]=0;
  }

  I[0]=N;


  for (int i=1; i<=N; i=i+1){

    //REAL mu=rand();

    //cout<<"pas N"<<i<<", rand="<<((REAL)rand())/RAND_MAX<<endl;

    REAL L=-log(((REAL)rand())/RAND_MAX)/st;
    X[i]=mu*L;
      //cout<<X[i]<<endl;


    //cout<<"pas N"<<i<<", X="<<X[i]<<endl;

    //pour avoir Q app exp(-st/nu)
    for (int j=1; j<=K; j++){
      if(X[i]>= (((REAL)j)/((REAL)(K))) ){
          //cout<<"itération"<<j<<endl;
          //cout<<X[i]<<"comparé"<<(((REAL)j)/((REAL)K))<<endl;
        I[j]=I[j]+1;
          //cout<<I[j]<<endl;
      }else break;
    }
  }

  for (int i=0; i<=K; i++){
    I[i]=I[i]/((REAL)(N*mu));
  }

  return I;
}

