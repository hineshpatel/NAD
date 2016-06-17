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
    for (auto step = 0; timeAcc < timeLimit; ++step, timeAcc += timeInc) {
        if ((activeBond.size()) || (np.position.z < (radius + bondCutoff.bondLMax))) {
            breakageCheck(); // assess bond breakage
            formationCheck(); // assess bond formation
            frepulsion = Frepulsion(np.position.z); // calculate repulsion force from substrate to nanoparticle
        }
        else frepulsion = coord{0, 0, 0};
        // don't have to check bond and calculate repulsion force
        fbond = Fbond(activeBond, bonds, receptors, ligands); // assess bond forces/torques on the nanoparticle
        fshear = Fshear(np.position.z); // assess shear force/torque on the nanoparticle
        acceleration(fbond, fshear, frepulsion); // calculate accelerations
        translation(np.velocity, np.position, np.acc); // translate nanoparticle
        rotateLig (ligands, rotate(np.rot_velocity, np.rot_acc), np.position); // rotate nanoparticle

        if (!(step%CHECKER)) {
            if (dist(np.lastPairPos,np.position)>bondL)getAvailRec(); // get all available receptors
            getAvailLig(); // get all available ligands
            if (ifDetach()) break;
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

//int main() {
//    setup();
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
//    set<int> activeBond = {0, 1, 2};
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
//    std::pair<coord, coord> bond = Fbond(activeBond, bonds, receptors, ligands);
//    std::pair<coord, coord> shear = Fshear(np.position.z);
//    acceleration(bond, shear, Frepulsion(np.position.z));
//    cout << Frepulsion(np.position.z) << "\t" << np.acc << "\t" << np.rot_acc << endl;
//
//}

/**
 * Note that
 * 1) the term "receptor" stands for adhesion molecule on the substrate;
 * 2) the term "ligand" stands for adhesion molecule on the nanoparticle;
 * 3) all counts start from 0;
 *
 **/