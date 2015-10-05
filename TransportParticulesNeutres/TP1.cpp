#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


using namespace std;

float Fuite(int N, float st, float x){

  float X[N];
  float mu=0.5;
  float Q=0;
  float u=0;

  srand (time(NULL));

  for (int i=1; i<=N; i=i+1){

      //float nu=rand();

      //cout<<"pas N"<<i<<", rand="<<((float)rand())/RAND_MAX<<endl;

      float L=-log(((float)rand())/RAND_MAX)/st;
      X[i]=mu*L;

      //cout<<"pas N"<<i<<", X="<<X[i]<<endl;

      //pour avoir Q app exp(-st/nu)
      if(X[i]>x){
        Q=Q+1;
      }
      //cout<<"pas N"<<i<<", Q="<<Q<<endl;
  }

  Q=Q/N;
  u=(Q/mu);

  return u;
}

float Fuite(int N, float st, float x, float mu){

  float X[N];

  float Q=0;
  float u=0;

  srand (time(NULL));

  for (int i=1; i<=N; i=i+1){

      //float nu=rand();

      //cout<<"pas N"<<i<<", rand="<<((float)rand())/RAND_MAX<<endl;

      float L=-log(((float)rand())/RAND_MAX)/st;
      X[i]=mu*L;

      //cout<<"pas N"<<i<<", X="<<X[i]<<endl;

      //pour avoir Q app exp(-st/nu)
      if(X[i]>x){
        Q=Q+1;
      }
      //cout<<"pas N"<<i<<", Q="<<Q<<endl;
  }

  Q=Q/N;
  u=(Q/mu);

  return u;
}
