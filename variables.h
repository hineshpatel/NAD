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

#define PI 3.14159

extern std::vector<receptor> receptors;
extern std::vector<ligand> ligands;
extern struct np np;
extern std::vector<bond> bonds;
extern std::set<int> activeBonds;
extern double timeAcc;
extern std::vector<int> availRec;
extern std::vector<int> availLig;
extern sfmt_t sfmt;


const double sigma_ts = _sigma; // (N/m)
const double volume = (4.0 / 3.0) * PI * _radius * _radius * _radius; // volume of NP (nm^3)
const double mass = _NPdensity * volume * 1.0e9; // mass of NP (ng)
const double rot_inertia = (2.0 / 5.0) * mass * _radius * _radius; // rotation inertia of NP (ng*nm^2)
const double beta = (3.0 * PI * _viscosity * (_radius * 2.0 * 0.001))
                    / (mass*1e-12); // see {English, Biophys.J, 2014} (s^-1)
const double beta_rot = PI * _viscosity * (2.0 * _radius * 0.001) * (2.0 * _radius) * (2.0 * _radius)
                        / rot_inertia * 1e12; // see {English, Biophys.J, 2014} (s^-1)
const double _substrate = 10 * _radius; // half diameter of the substrate (nm)
const int ligandNum = (int)(_ligandDens * 4 * PI * _radius * _radius / 1e6);
const int receptorNum = _receptorDens * _substrate * _substrate / 250000;
const unsigned long long steps = _timeLimit / _timeInc;
const cutoff bondCutoff{_bondEL, _bondCutRatio};



#endif //NANOAD_VARIABLES_H
