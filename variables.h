//
// Created by Mingqiu Wang on 5/7/16.
//


#ifndef NANOAD_VARIABLES_H
#define NANOAD_VARIABLES_H


#include <math.h>
#include <vector>
#include <map>
#include <set>
#include "dataStructures.h"
#include "parameters.h"
#include "SFMT.h"

//TODO: rearrange the variables
#define PI 3.14159

extern std::vector<receptor> receptors;
extern std::vector<ligand> ligands;
extern struct np np;
extern std::vector<bond> bonds;
extern std::set<int> activeBond;
extern double timeAcc;
extern std::vector<int> availRec;
extern std::vector<int> availLig;
extern cutoff bondCutoff;




const double sigma_ts = sigma; // (N/m)
const double volume = (4.0 / 3.0) * PI * radius * radius * radius; // volume of NP (nm^3)
const double mass = density * volume * 1.0e9; // mass of NP (ng)
const double rot_inertia = (2.0 / 5.0) * mass * radius * radius; // rotation inertia of NP (ng*nm^2)
const double beta = (3.0* PI * viscosity * (radius*2.0*0.001))
                    / (mass*1e-12); // see {English, Biophys.J, 2014} (s^-1)
const double beta_rot = PI*viscosity*(2.0*radius*0.001)*(2.0*radius)*(2.0*radius)
                        / rot_inertia*1e12; // see {English, Biophys.J, 2014} (s^-1)
const double substrate = 10*radius; // half diameter of the substrate (nm)
const int ligandNum = (int)(ligandDens * 4 * PI*radius*radius / 1e6);
const int receptorNum = receptorDens * substrate * substrate / 250000;
const unsigned long long steps = timeLimit/timeInc;



/******************* Random Number Generator **********************/
extern sfmt_t sfmt;


#endif //NANOAD_VARIABLES_H
