/**
 * Nano Adhesive Dynamics (NAD) simulation
 * A simulation framework for modeling adhesion/movement of adhesive nanoparticle in the flow
 *
 * Please modify "parameters.h" to customize simulation condition
 * Please refer to <> for simulation details
 *
 * @author: Mingqiu Wang (mingqiuw at uci.edu)
 * @date: 6/28/2016
 * @date: 5/04/2017 add attachment simulation.
 */

#include "declaration.h"

using namespace std;

int main() {


//    if (RESUME) { if (!resume()) exit(2); } // disable resume function in attachment sim
    ini();

// starts integrating Langevin equation
    coord frepulsion;
    pair<coord, coord> fbond, fshear;
    for (unsigned long long step = 0; timeAcc < _timeLimit; ++step, timeAcc += _timeInc) {

        if ((activeBonds.size()) || (np.position.z < (_radius + bondCutoff.bondLMax))) {
//            breakageCheck(activeBonds, bonds, ligands, receptors); // assess bond breakage
            if (formationCheck(availLig, availRec, activeBonds, ligands, receptors))
                attachedNP++;
                renewNP();
                // assess bond formation
            frepulsion = Frepulsion(np.position.z); // calculate repulsion force from substrate to nanoparticle
        }
        else frepulsion = coord{0, 0, 0}; // don't have to check bond and calculate repulsion force if
        // no bond and particle is higher than a certain height, which is not able to form bond
//        fbond = Fbond(activeBonds, bonds, receptors, ligands); // assess bond forces/torques on the nanoparticle
        fshear = Fshear(np.position.z); // assess shear force/torque on the nanoparticle
        acceleration(fbond, fshear, frepulsion); // calculate accelerations
        translation(np.velocity, np.position, np.acc); // translate nanoparticle
        if ((activeBonds.size()) || (np.position.z < (_radius + bondCutoff.bondLMax)))
            rotateLig (ligands, rotate(np.rot_velocity, np.rot_acc), np.position); // rotate nanoparticle

        if (!(step%CHECKER)) {
            if (np.position.z>=_boxHeight)
                renewNP();
            if (!inCell(np.position))
                putNPBack(np.position);

//            if (ifDetach(np.position)) break;
            checkDisplace(step);
        }
    }

// final recording
    writeEndTime();
//    writeBond();
    return 0;
}


/**
 * Note that
 * 1) the term "receptor" stands for adhesion molecule on the substrate;
 * 2) the term "ligand" stands for adhesion molecule on the nanoparticle;
 * 3) all counts start from 0;
 **/