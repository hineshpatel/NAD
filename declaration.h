//
// Created by Mingqiu Wang on 5/7/16.
//


#ifndef NANOAD_DECLARATION_H
#define NANOAD_DECLARATION_H

#include <stdio.h>
#include <stdlib.h>
#include <utility>
#include <iostream>
#include "outputs.h"
#include "parameters.h"
#include "SFMT.h"
#include "dataStructures.h"
#include "variables.h"

//TODO: rearrange the functions

void setup(); // pre-setup
double dist(double, double, double, double, double, double); // get distance between two points
double dist(const coord &coord1, const coord &coord2); // get distance between two points
std::pair<double, double> distribute
        (double upper_dis, double lower_dis, const std::pair<double, double>& i);
// distribute a point close (closer than upper_dis && further than lower_dis) to a given point i
void ini_receptor_double_cluster(); // initialize receptors
void ini_receptor_monomer(); // initialize receptors
void ini_receptor_single_cluster(); // initialize receptors
void ini_ligand(); // initialize ligands
void ini_binding(); // initialize binding
void ini_np(); // initialize nanoparticle
coord angle_trans(double ph, double tht, double radius); // transfer spherical to Cartesian
void getAvailRec();
void getAvailLig();
void breakageCheck();
void formationCheck();
coord Frepulsion(const double & npheight);
std::pair<coord, coord> Fshear(const double & npheight);
void acceleration(std::pair<coord, coord> Fbond, std::pair<coord, coord> Fshear, coord Frepulsion);
double gasdev(long idum);
std::pair<coord, coord> Fbond(const std::set<int> & activeBond, const std::vector<bond> & bonds,
                              const std::vector<receptor> & receptors, const std::vector<ligand> & ligands);
void translate();
void rotate();
int formBond(int lig, int rec);
void reporting();
bool ifDetach();
void writeBond();
void writeEndTime();
void writeLoc();
void writeResume();

#endif //NANOAD_DECLARATION_H
