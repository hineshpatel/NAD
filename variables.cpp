//
// Created by Mingqiu Wang on 5/8/16.
//
#include "variables.h"
#include "dataStructures.h"

std::vector<receptor> receptors;
std::vector<ligand> ligands;
sfmt_t sfmt;
struct np np;
double timeAcc;
std::vector<bond> bonds;
std::set<int> activeBond;
std::vector<int> availRec;
std::vector<int> availLig;
cutoff bondCutoff;
