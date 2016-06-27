/**
 * Nano Adhesive Dynamics (NAD) simulation
 * A simulation framework for modeling adhesion/movement of adhesive nanoparticle in the flow
 *
 * Please modify "parameters.h" to customize simulation condition
 * Please refer to <> for simulation details
 *
 * @author: Mingqiu Wang (mingqiuw at uci.edu)
 * @date: 6/23/2016
 *
 */

#include "declaration.h"

using namespace std;

int main() {

    if (RESUME) {
        if (!resume())
            exit(2);
    }
    else ini();

// get available adhesion molecules
    getAvailRec(availRec, np);
    getAvailLig(availLig, ligands);

// starts integrating Langevin equation
    coord frepulsion;
    pair<coord, coord> fbond, fshear;
    for (unsigned long long step = 0; timeAcc < timeLimit; ++step, timeAcc += _timeInc) {
        if ((activeBonds.size()) || (np.position.z < (_radius + bondCutoff.bondLMax))) {
            breakageCheck_linker(activeBonds, bonds, ligands, receptors); // assess bond breakage
            formationCheck(availLig, availRec, activeBonds, ligands, receptors); // assess bond formation
            frepulsion = Frepulsion(np.position.z); // calculate repulsion force from substrate to nanoparticle
        }
        else frepulsion = coord{0, 0, 0}; // don't have to check bond and calculate repulsion force if
        // no bond and particle is higher than a certain height, which is not able to form bond
        fbond = Fbond(activeBonds, bonds, receptors, ligands); // assess bond forces/torques on the nanoparticle
        fshear = Fshear(np.position.z); // assess shear force/torque on the nanoparticle
        acceleration(fbond, fshear, frepulsion); // calculate accelerations
        translation(np.velocity, np.position, np.acc); // translate nanoparticle
        rotateLig (ligands, rotate(np.rot_velocity, np.rot_acc), np.position); // rotate nanoparticle
//        printf("%e\t%lf\t%lf\t%lf\t%d\t%lf\t%lf\t%lf\n", timeAcc, np.position.x, np.position.y,
//               np.position.z, activeBonds.size(), np.acc.x, np.rot_acc.x, ligands.at(0).position.x);
//        if (step == 3400) {
//            exit(0);
//        }
        if (!(step%CHECKER)) {
            if (checkDisplace(step)) break;
        }

    }

// Step 8: final record
    writeEndTime();
    writeBond();
    return 0;
}

//int main() {
//    setRNG();
//    np.velocity = {28533362.197028,23205777.429526,-39590677.153833};
//    np.position = {1.553678, -45.818843, 124.235014};
////    np.acc = {3e9, 5e9, 9e7};
////    translation(np.velocity, np.position, np.acc);
////    cout << np.velocity << endl;
////    cout << np.position << endl;
//    np.rot_velocity = {478586.526767,-136219.076178,603880.867970};
//    np.rot_acc = {3e9, 5e9, 9e7};
//    ligand ligand1;
//    ligand1.updatePA({58.055188,-49.81867,35.823587}, np.position);
//    ligands.push_back(ligand1);
//    cout << "lig\t" << ligand1.position_origin << endl;
//    ligand1.updatePA({77.100255,15.987777,85.535451}, np.position);
//    ligands.push_back(ligand1);
//    cout << "lig\t" << ligand1.position_origin << endl;
//    ligand1.updatePA({-9.222984,-149.624129,135.781703}, np.position);
//    ligands.push_back(ligand1);
//    cout << "lig\t" << ligand1.position_origin << endl;
//    vector<coord> matrix = rotate(np.rot_velocity, np.rot_acc);
//    for (auto a : matrix)
//        cout << a << endl;
//    rotateLig (ligands, matrix, np.position);
//    for (auto a : ligands)
//        cout << a.position << endl;
//
//
//
//    }

//int main() {
//    set<int> activeBonds = {0, 1, 2};
//    vector<bond> bonds;
//    bond bond1;
//    bond1.ligand = 0;
//    bond1.receptor = 0;
//    bonds.push_back(bond1);
//    bond1.ligand = 1;
//    bond1.receptor = 1;
//    bonds.push_back(bond1);
//    bond1.ligand = 2;
//    bond1.receptor = 2;
//    bonds.push_back(bond1);
//    vector<receptor> receptors;
//    receptor receptor1;
//    receptor1.position = {0, -45, 0};
//    receptors.push_back(receptor1);
//    receptor1.position = {3, -50, 0};
//    receptors.push_back(receptor1);
//    receptor1.position = {-0.5, -40, 0};
//    receptors.push_back(receptor1);
//    np.position = {1.553678,-45.818843,124.235014};
//    vector<ligand> ligands;
//    ligand ligand1;
//    ligand1.updatePA({58.055188,-49.81867,35.823587}, np.position);
////    cout << ligand1.position_origin << endl;
//    ligands.push_back(ligand1);
//    ligand1.updatePA({77.100255,15.987777,85.535451}, np.position);
////    cout << ligand1.position_origin << endl;
//    ligands.push_back(ligand1);
//    ligand1.updatePA({-9.222984,-149.624129,135.781703}, np.position);
////    cout << ligand1.position_origin << endl;
//    ligands.push_back(ligand1);
//    std::pair<coord, coord> bond = Fbond(activeBonds, bonds, receptors, ligands);
//    std::pair<coord, coord> shear = Fshear(np.position.z);
//    acceleration(bond, shear, Frepulsion(np.position.z));
//    cout << Frepulsion(np.position.z) << "\t" << np.acc << "\t" << np.rot_acc << endl;
//
//}

//int main() {
//
//    FILE *outfile;
//    if ((outfile = fopen("test.txt", "a")) == NULL){ printf("\nerror on open test!"); exit(0); }
//
//    double length;
//    setRNG();
//    for (long x = 0; x<1e6; ++x) {
//        length = sfmt_genrand_res53(&sfmt)*5+38.6;
//        if (!ifBreak(length))
//            fprintf(outfile, "%lf\n", length);
//    }
//    fclose(outfile);
//
//}

/**
 * Note that
 * 1) the term "receptor" stands for adhesion molecule on the substrate;
 * 2) the term "ligand" stands for adhesion molecule on the nanoparticle;
 * 3) all counts start from 0;
 *
 **/