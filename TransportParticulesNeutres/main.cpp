#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include "TP1.h"
#include "SolveurMC.h"
#include "SolveurD.h"

using namespace std;

int main()
{
  char partie[]="Deterministe";//"Deterministe" ou "MonteCarlo"
  char materiau[]="Homogene"; //"Homogene" ou "Inhomogene" ou "HomogeneDiff"
  char source[]="Uniforme"; //"Dirac" ou "Uniforme" pour le solveurMCHomogene
  char methode[]="SN";
  REAL mu=0.5;


  //////////////////////////////////////////////////
 ////////////  Partie Monte Carlo ////////////////////
 ///////////////////////////////////////////////////

 if(strcmp(partie,"MonteCarlo")==0){

    //parametres du solveur
    int N=1000000;//nombre de particules
    int K=100;//taille de la discrétisation de [0,1]

    //Initialistaion des constantes
    //On défini les caractéristiques du problème

    //definition du materieu Inhomogene
    std::vector<REAL> abst (3);
    std::vector<REAL> st (3);
    abst[0]=0;
    st[0]=1;
    abst[1]=0.3;
    st[1]=3;
    abst[2]=0.7;
    st[2]=1;
    int sizest=st.size();

    //déclaration des fichiers texte de sorti
    std::ofstream solveur ("solveurMC.dat");
    std::ofstream exact ("exactMC.dat");


    std::vector<REAL> absm (K+1);
    std::vector<REAL> E (K+1);
    std::vector<REAL> Im (K+1);

    //initialisation de l'axe des abscisses
    for (int i=0 ; i<absm.size() ; ++i) {
          absm[i] = i/((REAL)K);
    }

    //Appel aux solveurs
      //Homogène
    if (strcmp(materiau,"Homogene")==0){
        std::vector<REAL> Im=SolveurMCHomogene(N, K, 1, mu, source);
        //écriture dans un fichier texte
        for (int i=0 ; i<absm.size() ; ++i) {
          solveur << absm[i] << " " << Im[i] << std::endl;
        }
    }
      //Inhomogène
    if (strcmp(materiau,"Inhomogene")==0){
        std::vector<REAL> Im=SolveurMCNonHomogene (N, K, mu, abst, st, sizest);
        //écriture dans un fichier texte
        for (int i=0 ; i<absm.size() ; ++i) {
          solveur << absm[i] << " " << Im[i] << std::endl;
        }
    }

    if (strcmp(materiau,"HomogeneDiff")==0){
        std::vector<REAL> Im=SolveurMCHomogeneDiff (N, K);
        //écriture dans un fichier texte
        for (int i=0 ; i<absm.size() ; ++i) {
          solveur << absm[i] << " " << Im[i] << std::endl;
        }
    }

    //Calcul de la solution exacte
    if (strcmp(materiau,"Homogene")==0){
      if (strcmp(source,"Dirac")==0){
        for (int i=0 ; i<absm.size() ; ++i) {
          //solution exacte pour source dirac 0
          E[i]=(1/mu)*exp(-1*absm[i]/mu);
        }
      }
      if (strcmp(source,"Uniforme")==0){
        if(mu>=0){
          for (int i=0 ; i<absm.size() ; ++i) {
            //solution exacte pour source uniforme, mu positif
            E[i]=(1/1)*(1-exp(-1*absm[i]/mu));
          }
        }
        if(mu<0){
          for (int i=0 ; i<absm.size() ; ++i) {
            //solution exacte pour source uniforme, mu négatif
            E[i]=(1/1)*(1-exp(-1*(absm[i]-1)/mu));
          }
        }
      }
    }


    if (strcmp(materiau,"Inhomogene")==0){
      for (int i=0 ; i<absm.size() ; ++i) {
         //solution exacte pour source uniforme, mu négatif
         if ( absm[i] <= 0.3 ){
           E[i]=(1/mu)*exp(-absm[i]/mu);
         }
         if ( (0.3<absm[i])&&(absm[i]<=0.7) ){
           E[i]=(1/mu)*exp(-1*(3*absm[i]-0.6)/mu);
         }
         if ( 0.7<absm[i] ){
           E[i]=(1/mu)*exp(-1*(absm[i]+0.8)/mu);
         }
      }
    }


    if (strcmp(materiau,"HomogeneDiff")==0){
      for (int i=0 ; i<absm.size() ; ++i) {
         //On ne connait pas la solution
          E[i]=0;
      }
    }




    //Ecritures de la solution exacte résultats
    for (int i=0 ; i<absm.size() ; ++i) {
      exact << absm[i] << " " << E[i]<< std::endl;
    }
    //cout<<"u="<<u<<endl;
    //cout<<"u réél="<<2*exp(-2)<<endl;
 }


 ////////////////////////////////////////////////////////
//////////// ~~~~~> Partie Deterministe <~~~~~ ////////////
 ///////////////////////////////////////////////////////

  if ( (strcmp(partie,"Deterministe")==0) && (strcmp(methode,"SN")>0)  ){

    int K=10000;


    std::vector<REAL> sd (K+1);// <3   <3    <3
    std::vector<REAL> Qt (K+1);//   <3    <3
    std::vector<REAL> absm (K+1);
    REAL CL=1/mu;//1/mu, mu>0; 0 sinon
    if(strcmp(source,"Uniforme")==0){
      CL=0;
    }

    std::vector<REAL> I (K+1);

    for (int i=0; i<K+1; i++){
      I[i]=0;
      Qt[i]=0;
      if(strcmp(source,"Uniforme")==0){
        Qt[i]=1;
      }
      absm[i] = i/((REAL)K);
      sd[i]=1;
      if (strcmp(materiau,"Inhomogene")==0){
        if ( (0.3<absm[i])&&(absm[i]<=0.7) ){
          sd[i]=3;
        }
      }
    }


    std::ofstream solveur ("solveurD.dat");

    I=SchemaDiamant(K,mu,sd,Qt,CL);

    for(int i=0;i<K+1;i++){
      solveur<<absm[i]<<" "<<I[i]<<std::endl;
    }
  }

  if ( (strcmp(partie,"Deterministe")==0) && (strcmp(methode,"SN")==0)  ){
    int K=10000;
    REAL CL;
    REAL epsilon=0.01;
    std::vector<REAL> S (K);
    std::vector<REAL> sd (K);
    std::vector<REAL> absm (K);

    int i=0;

    for (i=0; i<K; i++){
      sd[i]=1;
      absm[i]=(1/2*K)+((REAL)i)/K;
    }

    if (strcmp(source,"Uniforme")==0){
      CL=0;
      for (i=0; i<K+1; i++){
        S[i]=1;
      }
    }
    if (strcmp(source,"Dirac")==0){
      CL=0;
      for (i=0; i<K+1; i++){
        S[i]=0;
      }
      S[0]=1;
    }

    std::vector<REAL> I (K);
    std::ofstream solveur ("solveurD.dat");

    I=MethodeSN( K, epsilon, CL, sd, S);

    for(int i=0;i<K;i++){
      solveur<<absm[i]<<" "<<I[i]<<std::endl;
    }
  }


}
