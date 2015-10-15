#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include "SolveurMC.h"


using namespace std;

std::vector<REAL> SolveurMCHomogene(int N, int K, REAL st, REAL mu, char* source){


  //Initialisation des constantes

  std::vector<REAL> X (N);
  std::vector<REAL> I(K+1);
  REAL X0=0;
  int j;

  srand (time(NULL));


  for (int i=0; i<K; i=i+1){
    I[i]=0;
  }

  for (int i=1; i<=N; i=i+1){


    REAL L=-log(((REAL)rand())/RAND_MAX)/st;

    if (strcmp(source,"Dirac")==0){
      X0=0;
      I[0]=N;
    }
    else if (strcmp(source,"Uniforme")==0){
      X0= ((REAL)rand())/RAND_MAX;
    }

    X[i]=X0+mu*L;


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

std::vector<REAL> SolveurMCNonHomogene(int N, int K,REAL mu, vector<REAL> abst, vector<REAL> st, int sizest){


  //Initialisation des constantes

  std::vector<REAL> X (N);
  std::vector<REAL> I(K+1);
  REAL X0=0;
  int j;
  int i=0; //marque la particule en cours de traitement
  int k;//marque l'intervalle abst ou se trouve X[i]

  //Initialisation à 0 des vecteurs de resultats
  srand (time(NULL));
  for (int i=0; i<N; i=i+1){
    X[i]=0;
  }
  for (int i=0; i<K+1; i=i+1){
    I[i]=0;
  }


  REAL stmax=0;
  for(int i=0;i<sizest;i++){
    if(stmax<st[i]){
      stmax=st[i];
    }
  }


  //Traitement des particules
  i=0;
  while(i<N){

    X0=X[i];

    REAL L=-log(((REAL)rand())/RAND_MAX)/stmax;

    X[i]=X[i]+mu*L;


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
        //cout<<"J'incrémente la frontière"<<(((REAL)j)/((REAL)(K)))<<endl;
        }
      }
    }

    k=0;
    j=0;

    for (int j=0; j<sizest;j++){
      if(X[i]>=abst[j]){
        k=j;
      }
    }
    //choc fictif ?
    if(  (((REAL)rand())/RAND_MAX)  < (st[k]/stmax) ){
      i=i+1;//absorbtion
    }
  }

  for (int i=0; i<=K; i++){
    I[i]=fabs(I[i]/((REAL)(N*mu)));
    //cout<<I[i]<<endl;
  }

  return I;
}

std::vector<REAL> SolveurMCHomogeneDiff(int N, int K){


  //Initialisation des constantes

  std::vector<REAL> X (N);
  std::vector<REAL> I(K+1);
  REAL X0=0;
  int j;
  int i=0; //marque la particule en cours de traitement
  REAL st=1;
  REAL ss=1;

  //Initialisation à 0 des vecteurs de resultats
  srand (time(NULL));
  for (int i=0; i<N; i=i+1){
    X[i]=((REAL)rand())/RAND_MAX;//initilisation aléatoire des positions entre 0 et 1 (source uniforme)
  }
  for (int i=0; i<K+1; i=i+1){
    I[i]=0;
  }


  //Traitement des particules
  i=0;
  while(i<N){

    X0=X[i];

    REAL L=-log(((REAL)rand())/RAND_MAX)/st;
    REAL mu=((REAL)rand())/RAND_MAX;

    X[i]=X[i]+mu*L;


    for (int j=0; j<K+1; j++){
        if( (X0<=(((REAL)j)/((REAL)(K))))&&((((REAL)j)/((REAL)(K)))<=X[i]) ){
          I[j]=I[j]+(1/mu);
        }
    }
    j=0;

    //Diffusion ?
    if(  (((REAL)rand())/RAND_MAX)  <=  (st/(st+ss))   ){
      i=i+1;//on passe a la particule suivante
    }
  }

  for (int i=0; i<=K; i++){
    I[i]=fabs(I[i]/((REAL)(N)));
    //cout<<I[i]<<endl;
  }

  return I;
}
