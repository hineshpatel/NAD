////
//// Created by Mingqiu Wang on 6/29/16.
////
//
///**
// * Nano Adhesive Dynamics (NAD) simulation
// * A simulation framework for modeling adhesion/movement of adhesive nanoparticle in the flow
// *
// * For modeling Linker/Ab/ICAM-1
// *
// * Please modify "parameters.h" to customize simulation condition
// * Please refer to <> for simulation details
// *
// * @author: Mingqiu Wang (mingqiuw at uci.edu)
// * @date: 6/29/2016
// *
// */
//
//#include "declaration.h"
//
//using namespace std;
//
//int main() {
//
//    cout << beta << "\t" << mass <<  endl;
//    cout << beta_rot << "\t" << rot_inertia <<  endl;
//
//    if (RESUME) { if (!resume()) exit(2); }
//    else ini();
//
//// starts integrating Langevin equation
//    coord frepulsion;
//    pair<coord, coord> fbond, fshear;
//    for (unsigned long long step = 0; timeAcc < _timeLimit; ++step, timeAcc += _timeInc) {
//
//        if ((activeBonds.size()) || (np.position.z < (_radius + bondCutoff.bondLMax))) {
//            breakageCheck_linker(activeBonds, bonds, ligands, receptors); // assess bond breakage
//            formationCheck(availLig, availRec, activeBonds, ligands, receptors); // assess bond formation
//            frepulsion = Frepulsion(np.position.z); // calculate repulsion force from substrate to nanoparticle
//        }
//        else frepulsion = coord{0, 0, 0}; // don't have to check bond and calculate repulsion force if
//        // no bond and particle is higher than a certain height, which is not able to form bond
//        fbond = Fbond_linker(activeBonds, bonds, receptors, ligands); // assess bond forces/torques on the nanoparticle
//        fshear = Fshear(np.position.z); // assess shear force/torque on the nanoparticle
//        acceleration(fbond, fshear, frepulsion); // calculate accelerations
//        translation(np.velocity, np.position, np.acc); // translate nanoparticle
//        rotateLig (ligands, rotate(np.rot_velocity, np.rot_acc), np.position); // rotate nanoparticle
//        if (!(step%CHECKER)) {
//            if (ifDetach(np.position)) break;
//            checkDisplace(step);
//        }
//    }
//
//// final recording
//    writeEndTime();
//    writeBond();
//    return 0;
//}
//
//
///**
// * Note that
// * 1) the term "receptor" stands for adhesion molecule on the substrate;
// * 2) the term "ligand" stands for adhesion molecule on the nanoparticle;
// * 3) all counts start from 0;
// **/