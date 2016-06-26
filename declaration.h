//
// Created by Mingqiu Wang on 5/7/16.
//


#ifndef NANOAD_DECLARATION_H
#define NANOAD_DECLARATION_H

#include <stdio.h>
#include <stdlib.h>
#include <utility>
#include <iostream>
#include <queue>
#include <fstream>

#include "outputs.h"
#include "parameters.h"
#include "SFMT.h"
#include "dataStructures.h"
#include "variables.h"

//TODO: rearrange the functions

void setRNG(sfmt_t & sfmt);
double dist(double, double, double, double, double, double); // get distance between two points
double dist(const coord &coord1, const coord &coord2); // get distance between two points
std::pair<double, double> distribute
        (double upper_dis, double lower_dis, const std::pair<double, double>& i);
// distribute a point close (closer than upper_dis && further than lower_dis) to a given point i
void ini_receptors_doubleCluster(std::vector<receptor> &receptors); // initialize receptors
void ini_receptor_monomer(std::vector<receptor> & receptors); // initialize receptors
void ini_receptor_cluster(std::vector<receptor> &receptors); // initialize receptors
void ini_ligand(std::vector<ligand> & ligands); // initialize ligands
void ini_binding(std::vector<receptor> & receptors, std::vector<ligand> & ligands,
                 std::set<int> & activeBonds, std::vector<bond> & bonds, const struct np & np); // initialize binding
void ini_np(struct np & np); // initialize nanoparticle
coord angle_trans(double ph, double tht, double radius); // transfer spherical to Cartesian
void getAvailRec(std::vector<int> & availRec, struct np & np);
void getAvailLig(std::vector<int> & availLig, const std::vector<ligand> & ligands);
void breakageCheck(std::set<int> &activeBonds, std::vector<bond> & bonds,
                   std::vector<ligand> & ligands, std::vector<receptor> & receptors);
void formationCheck(const std::vector<int> & availLig, const std::vector<int> & availRec,
                    std::set<int> &activeBonds, std::vector<ligand> & ligands,
                    std::vector<receptor> & receptors);
coord Frepulsion(const double & npheight);
std::pair<coord, coord> Fshear(const double & npheight);
void acceleration(std::pair<coord, coord> Fbond, std::pair<coord, coord> Fshear, coord Frepulsion);
double gasdev(long idum);
std::pair<coord, coord> Fbond(const std::set<int> & activeBond, const std::vector<bond> & bonds,
                              const std::vector<receptor> & receptors, const std::vector<ligand> & ligands);
void translation(coord &velocity, coord &position, const coord &acc);
std::vector<coord> rotate(coord & velocity, const coord & acc);
void rotateLig (std::vector<ligand> & ligands, const std::vector<coord> & rotationMatrix, const coord & NPposition);
int formBond(int lig, int rec, std::vector<ligand> & ligands, std::vector<receptor> & receptors,
             std::vector<bond> & bonds);
bool ifBreak(double bondLength);
bool checkDisplace(unsigned long long &step);
bool ifDetach(const coord & position);
void writeBond();
void writeEndTime();
void writeLoc();
void writeResume();
double distSQ(const coord &coord1, const coord &coord2);
bool ifForm(double delta);
double distance_P_S(const coord & a0, const coord & a1, const coord & a2);
double seg_seg_Dist (const coord & rec1, const coord & lig2, const coord & rec3, const coord & lig4);
bool ifCross (const std::set<int> & activeBond, const std::vector<receptor> & receptors,
              const std::vector<ligand> & ligands, const int checkLig, const int checkRec);
double dist(const std::pair<double, double> & coord1, const std::pair<double, double> & coord2);

double delta2_com(double);
double delta2_ext(double);
void breakageCheck_linker(std::set<int> &activeBonds, std::vector<bond> & bonds,
                          std::vector<ligand> & ligands, std::vector<receptor> & receptors);
bool resume();
void ini();
#endif //NANOAD_DECLARATION_H
