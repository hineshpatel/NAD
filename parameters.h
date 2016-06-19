//
// Created by Mingqiu Wang on 5/7/16.
//

#ifndef NANOAD_PARAMETERS_H
#define NANOAD_PARAMETERS_H

//TODO: rearrange the parameters

#include <math.h>

/* =======  Commonly used  ====================================================== */
const int receptorDens = 21; // ICAM-1 density (1/um2)
const int ligandDens = 410; // Ab density (1/um2)
const double bondL = 41.1; // bond length at equilibrium (nm)
const double timeLimit = 5; // (s)
const double sigma = 0.8; // (N/m)
const double gama = 0.274; // (nm)
#define LIG_CLU 1
#define REC_CLU 1
#define DOUBLE_CLU 1

/* =======  Particle properties  ====================================================== */
const double radius = 105; // radius of the NP (nm)
const double density = 1.05e-21; // density of polysterene (g/nm^3)


/* =======  Computational properties  ====================================================== */
const double timeInc = 1.0e-9; // step length (s)


/* =======  Physical properties  ====================================================== */
#define BROWNIAN 1
const double viscosity = 1e-9; // viscosity of water (kg/um/s)
const double thermal = 4.11; // scaled temperature Kb*T (nm^2*g/s^2)

/* =======  For receptor double clustering  ====================================================== */
const double clus_maxRecDis = 2; // max distance between 2 clustered receptors (nm)
const double clus_minRecDis = 0.5; // min distance between 2 clustered receptors (nm)
const double clus_gProDis = 4; // max distance from a receptor to a G protein (nm)
const int gProteinDens = 134/4 + 1; // density of G protein sites (1/um2)

/* =======  For receptor clustering  ====================================================== */

/* =======  For ligand clustering  ====================================================== */
const double clus_maxLigDis = 10; // max distance between 2 fab domain of an Ab (nm)
const double clus_minLigDis = 5; // min distance between 2 fab domain of an Ab (nm)

const double kr0 = 1.1e-4; // (s^-1)
const double kf0 = 1.6e5; // (s^-1)
const double shear_rate = 100; // (s^-1)
#define DEBUG 1
const double proThick = 5; // the combined thickness of the protein layers (nm)
const double compressibility = 3e-5; // the compressibility coefficient of surface proteins (nN);
const int detach_criteria = 5; // # times than radius

#define REPORTER 10000000  // all outputs saved every * steps
#define CHECKER 100
#define CROSSCHECK 1 // whether to perform bond cross check

const double bondDiameter = 2; // ICAM-1/Ab bond diameter (nm)
#endif //NANOAD_PARAMETERS_H
