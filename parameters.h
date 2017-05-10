//
// Created by Mingqiu Wang on 5/7/16.
//

#ifndef NANOAD_PARAMETERS_H
#define NANOAD_PARAMETERS_H

#include <math.h>
#include "dataStructures.h"

/* =======  Commonly used  ====================================================== */

const int _ligandDens = 3400; // Ab density (1/um2)
const int _receptorDens = 134; // ICAM-1 density (1/um2)
const double _bondEL = 41.1; // bond length at equilibrium (nm)
const double _timeLimit = 30; // (s)
const double _sigma = 0.001; // stiffness of ICAM-1/Ab bond (N/m)
const double _bondCutRatio = 0.5;
const double sigma_ts = 0.00001; // (N/m)
const double _gama = 11; // reactive compliance of ICAM-1/Ab bond (nm)
#define LIG_CLU 1 // if cluster two binding sites of Ab
#define REC_CLU 1 // if cluster two ICAM-1
#define DOUBLE_CLU 1 // if double clustering four ICAM-1


/* ======= Output ============================================================== */
#define WRITEBONDF 0
#define WRITEBONDL 0

/* =======  Particle properties  =============================================== */

const double _radius = 105; // _radius of the NP (nm)
const double _NPdensity = 1.05e-21; // density of polysterene (g/nm^3)


/* ======= Attachment ========================================================== */
const double _boxHeight = 300; // the height of box in which NP is moving (nm)
const double _boxLength= 8 * _radius; // the x/y radius of box in which NP is moving (nm)


/* =======  Computational properties  ========================================== */

const double _timeInc = 1.0e-9; // step length (s)
const int _detach_criteria = 5; // # times than _radius
#define DEBUG 1
#define RESUME 0
#define REPORTER 100000000  // all outputs saved every * steps
#define TRAJ 10000000  // np location saved every * steps
#define CHECKER 100


/* =======  Physical properties  =============================================== */

#define BROWNIAN // If with Brownian motion or not
const double _viscosity = 1e-9; // viscosity of water (kg/um/s)
const double _thermal = 4.11; // scaled temperature Kb*T (nm^2*g/s^2)
const double _kr0 = 1.1e-4; // (s^-1)
const double _kf0 = 1.6e5; // (s^-1)
const double _shear_rate = 100; // (s^-1)
const double _proThick = 5; // the combined thickness of the protein layers (nm)
const double _compressibility = 3e-5; // the _compressibility coefficient of surface proteins (nN);
#define CROSSCHECK 1 // whether to perform bond cross check
const double _bondDiameter = 2; // ICAM-1/Ab bond diameter (nm)


/* =======  For receptor double clustering  ==================================== */

const double _clus_gProDis = 4; // max distance from a receptor to a G protein (nm)
const int _gProteinDens = 134 / 4 + 1; // density of G protein sites (1/um2)


/* =======  For receptor clustering  =========================================== */

const double _clus_maxRecDis = 2; // max distance between 2 clustered receptors (nm)
const double _clus_minRecDis = 0.5; // min distance between 2 clustered receptors (nm)


/* =======  For ligand clustering  ============================================= */

const double clus_maxLigDis = 10; // max distance between 2 fab domain of an Ab (nm)
const double clus_minLigDis = 5; // min distance between 2 fab domain of an Ab (nm)


/* =======  For springy linker  ================================================ */
const double _linkerLen = 2.6; // (nm)
const double _maxDeltaCom = _linkerLen - 0.2; // (nm)
const double _maxDeltaExt = 13 - _linkerLen;
const double _sigma_linker = 1e-3; // (N/m)


#endif //NANOAD_PARAMETERS_H
