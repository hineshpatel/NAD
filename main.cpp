/**
 * Nano Adhesive Dynamics (NAD) simulation
 * A simulation framework for modeling adhesion/movement of adhesive nanoparticle in the flow.
 *
 * Please modify "parameters.h" to customize simulation configurations.
 * Please refer to <Wang, Mingqiu, et al. Langmuir 32.49 (2016): 13124-13136.>
 * for simulation details.
 *
 * @author: Mingqiu Wang (mingqiuw at uci.edu)
 * @date: 6/2/2017
 * @date: 11/2/2017 (merged detachment sim and attachment sim)
 */

#include "declaration.h"

using namespace std;

int att() {
    ini();
    // starts integrating Langevin equation
    coord frepulsion;
    pair<coord, coord> fbond, fshear;
    for (unsigned long long step = 0; timeAcc < _timeLimit; ++step, timeAcc += _timeInc) {

        if ((activeBonds.size()) || (np.position.z < (_radius + bondCutoff.bondLMax))) {
            bool formed;
            if (ORI)
                formed = formationCheckOri(availLig, availRec, activeBonds, ligands, receptors);
            else
                formed = formationCheck(availLig, availRec, activeBonds, ligands, receptors);

            if (formed) {
                writeResume(attachedNP++);
                writeBond(attachedNP);
                activeBonds.clear();
                renewNP();
            }
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

int detach() {
    if (RESUME) { if (!resume()) exit(2); }
    else ini();

    // starts integrating Langevin equation
    coord frepulsion;
    pair<coord, coord> fbond, fshear;
    for (unsigned long long step = 0; timeAcc < _timeLimit; ++step, timeAcc += _timeInc) {

        if ((activeBonds.size()) || (np.position.z < (_radius + bondCutoff.bondLMax))) {
            breakageCheck(activeBonds, bonds, ligands, receptors); // assess bond breakage
            if (ORI)
                formationCheckOri(availLig, availRec, activeBonds, ligands, receptors);
            else
                formationCheck(availLig, availRec, activeBonds, ligands, receptors);
            // assess bond formation
            frepulsion = Frepulsion(np.position.z); // calculate repulsion force from substrate to nanoparticle
        }
        else frepulsion = coord{0, 0, 0}; // don't have to check bond and calculate repulsion force if
        // no bond and particle is higher than a certain height, which is not able to form bond
        fbond = Fbond(activeBonds, bonds, receptors, ligands); // assess bond forces/torques on the nanoparticle
        fshear = Fshear(np.position.z); // assess shear force/torque on the nanoparticle
        acceleration(fbond, fshear, frepulsion); // calculate accelerations
        translation(np.velocity, np.position, np.acc); // translate nanoparticle
        rotateLig (ligands, rotate(np.rot_velocity, np.rot_acc), np.position); // rotate nanoparticle
        if (!(step%CHECKER)) {
            if (ifDetach(np.position)) break;
            checkDisplace(step);
        }
    }

    // final recording
    writeEndTime();
    writeBond();
    return 0;
}
int main() {
    if (ATT) {
        att();
    } else {
        detach();
    }
}


/**
 * Note that
 * 1) the term "receptor" stands for adhesion molecule on the substrate;
 * 2) the term "ligand" stands for adhesion molecule on the nanoparticle;
 * 3) all counts start from 0;
 **/