#ifndef SOLVEURMC_H
#define SOLVEURMC_H
#include<vector>

typedef double REAL;

std::vector<REAL> SolveurMCHomogene(int, int, REAL, REAL, char*);

std::vector<REAL> SolveurMCNonHomogene(int, int, REAL, std::vector<REAL>, std::vector<REAL>, int);

std::vector<REAL> SolveurMCHomogeneDiff(int, int);


#endif

