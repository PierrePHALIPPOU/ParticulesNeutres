#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include<vector>
#include "SolveurMC.h"


using namespace std;

std::vector<REAL> SolveurMC(int N, int K, REAL st, REAL mu, char* source){

  std::vector<REAL> X (N);
  std::vector<REAL> I(K+1);
  REAL Q=0;
  REAL u=0;
  REAL X0=0;
  int j;

  srand (time(NULL));

//initialisation de I

  for (int i=0; i<K; i=i+1){
    I[i]=0;
  }


  //for (int i=1; i<=1; i=i+1){
  for (int i=1; i<=N; i=i+1){


    REAL L=-log(((REAL)rand())/RAND_MAX)/st;

    if (strcmp(source,"Dirac")==0){
      X0=0;
      I[0]=N;
    }
    else if (strcmp(source,"Uniforme")==0){
      X0= ((REAL)rand())/RAND_MAX;
      //cout << "XO="<<X0<<endl;
    }

    X[i]=X0+mu*L;
    //cout<<"X="<<X[i]<<endl;


    for (int j=0; j<K+1; j++){
      if(mu>=0){
        if( (X0<=(((REAL)j)/((REAL)(K))))&&((((REAL)j)/((REAL)(K)))<=X[i]) ){
        I[j]=I[j]+1;
        //cout<<"J'incrémente la frontière"<<(((REAL)j)/((REAL)(K)))<<endl;
        }
      }
      if(mu<0){
        if( (X[i]<=(((REAL)j)/((REAL)(K))))&&((((REAL)j)/((REAL)(K)))<=X0) ){
        I[j]=I[j]+1;
        }
      }

    }
  }

  for (int i=0; i<=K; i++){
    I[i]=fabs(I[i]/((REAL)(N*mu)));
  }

  return I;
}

