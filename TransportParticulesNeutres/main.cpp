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
      //On défini les caractéristiques du problème
      char source[]="Uniforme"; //"Dirac" ou "Uniforme" pour le solveurMCHomogene
      char materiau[]="HomogeneDiff"; //"Homogene" ou "Inhomogene" ou "HomogeneDiff"


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
    std::ofstream solveur ("solveur.dat");
    std::ofstream exact ("exact.dat");

    //parametres du solveur
    int N=1000000;//nombre de particules
    int K=100;//taille de la discrétisation de [0,1]
    REAL mu=0.5;//mu (non utilisé dans solveurMCHomogeneDiff

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
