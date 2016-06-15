#include <iostream>
#include "declaration.h"

//TODO: restart
using namespace std;


int main() {

// Step 1: sets up the simulation
    setup();

// Step 2: adds receptors to substrate
    if (REC_CLU&&DOUBLE_CLU)
        ini_receptor_double_cluster(); // as double clustering
    else if (REC_CLU)
        ini_receptor_single_cluster(); // as single clustering
    else ini_receptor_monomer(); // as monomer

// Step 3: adds a nanoparticle to the system
    ini_np();

// Step 4: adds ligands on the nanoparticle
    ini_ligand();

// Step 5: sets up initial NP binding
    ini_binding();

// Step 6: get available adhesion molecules
    getAvailRec();
    getAvailLig();


// Step 7: starts integrating Langevin equation
    coord frepulsion;
    pair<coord, coord> fbond, fshear;
    for (auto step = 0; step < steps; ++step, timeAcc += timeInc) {
        if ((activeBond.size()) || (np.position.z < (radius + bondCutoff.bondLMax))) {
            breakageCheck(); // assess bond breakage
            formationCheck(); // assess bond formation
            frepulsion = Frepulsion(); // calculate repulsion force from substrate to nanoparticle
        }
        else frepulsion = coord{0, 0, 0};
        // don't have to check bond and calculate repulsion force
        fbond = Fbond(); // assess bond forces/torques on the nanoparticle
        fshear = Fshear(); // assess shear force/torque on the nanoparticle
        acceleration(fbond, fshear, frepulsion); // calculate accelerations
        translate(); // translate nanopartile
        rotate(); // rotate nanoparticle

        if (!(step%CHECKER)) {
            if (dist(np.lastPairPos,np.position)>bondL)getAvailRec(); // get all available receptors
            getAvailLig(); // get all available ligands
            if (ifDetach()) break;
            reporting();
        }
        if (!(step%REPORTER)) {
            writeBond();
            writeLoc();
            writeResume();
        }
    }

// Step 8: final record
    writeEndTime();



    return 0;
}


/**
 * Note that
 * 1) the term "receptor" stands for adhesion molecule on the substrate;
 * 2) the term "ligand" stands for adhesion molecule on the nanoparticle;
 * 3) all counts start from 0;
 *
 **/