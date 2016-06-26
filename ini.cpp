//
// Created by Mingqiu Wang on 5/7/16.
//

#include <fstream>
#include "declaration.h"
using namespace std;
/**
 * This function initializes receptors on the substrate under double clustering condition.
 *
 * @update: receptors: initialize all receptors
 *
 */
void ini_receptor_double_cluster(std::vector<receptor> & receptors) {

    int gProteinNum = _gProteinDens * _substrate * _substrate / 250000; // number of G protein sites
    std::vector<int> gProtein (gProteinNum + 1, 0); // how many receptors in each G protein site (max 4)
    std::vector<std::pair<double, double> > gProteinPos; // positions of G protein sites (nm)
    double x, y; srand(sfmt_genrand_res53(&sfmt));

    // distributes G protein sites on the substrate
    for (auto j = 0; j < gProteinNum; j++) {
        x = sfmt_genrand_res53(&sfmt) * 2 * _substrate - _substrate;
        y = sfmt_genrand_res53(&sfmt) * 2 * _substrate - _substrate;
        gProteinPos.push_back({x, y });
    }


    int i;
    double receptor_x, receptor_y;

    // distributes receptors
    for (auto j = 0; j < receptorNum; j++) {
        if (j % 2) {
            std::pair<double, double> a = distribute(clus_maxRecDis,
                                                     clus_minRecDis, {receptors.at(j - 1).position.x, receptors.at(j - 1).position.y });
            receptor_x = a.first;
            receptor_y = a.second;
        }
        else {
            do {
                i = rand() % gProteinNum + 1;
            } while (gProtein[i] == 4);
            std::pair<double, double> a = distribute(clus_gProDis, 0, gProteinPos[i]);
            receptor_x = a.first;
            receptor_y = a.second;
        }
        receptor receptor1;
        receptor1.position = {receptor_x, receptor_y, 0};
        gProtein[i]++;
        receptors.push_back(receptor1);

    }
}


/**
 * This function initializes receptors on the substrate under single clustering condition.
 *
 * @update: receptors: initialize all receptors
 *
 */
void ini_receptor_single_cluster(std::vector<receptor> & receptors) {

    double receptor_x, receptor_y;
    for (auto j = 0; j < receptorNum; j++) { // distribute receptors
        if (j % 2) {
            std::pair<double, double> a = distribute(clus_maxRecDis,
                                                     clus_minRecDis, {receptors.at(j - 1).position.x, receptors.at(j - 1).position.y });
            receptor_x = a.first;
            receptor_y = a.second;
        }
        else {
            receptor_x = sfmt_genrand_res53(&sfmt) * 2 * _substrate - _substrate;
            receptor_y = sfmt_genrand_res53(&sfmt) * 2 * _substrate - _substrate;
        }
        receptor receptor1;
        receptor1.position = {receptor_x, receptor_y, 0};
        receptors.push_back(receptor1);
    }
}


/**
 * This function initializes receptors on the substrate under monomer condition.
 *
 * @update: receptors: initialize all receptors
 *
 *
 */
void ini_receptor_monomer(std::vector<receptor> & receptors) {

    double receptor_x, receptor_y;
    for (auto j = 0; j < receptorNum; j++) { // distribute receptors

        receptor_x = sfmt_genrand_res53(&sfmt) * 2 * _substrate - _substrate;
        receptor_y = sfmt_genrand_res53(&sfmt) * 2 * _substrate - _substrate;

        receptor receptor1;
        receptor1.position = {receptor_x, receptor_y, 0};
        receptors.push_back(receptor1);
    }

}


/**
 * This function initializes a nanoparticle.
 *
 * @update: np
 *
 */
void ini_np(struct np & np) {

    np.position.x = 0;
    np.position.y = 0;
    np.position.z = _radius + _bondEL;
    np.lastPairPos = np.position;

}


/**
 * This function initializes ligandss on the NP.
 *
 * @update: ligands
 *
 */
void ini_ligand(std::vector<ligand> & ligands) {

    double ph, tht, dis; coord lig_ori;
    for (auto j = 0; j < ligandNum; j++) {
        if (LIG_CLU) {
            if (j % 2) {
                do {
                    ph = 2.0 * PI * sfmt_genrand_res53(&sfmt); // [0, 2PI]
                    tht = acos(2 * sfmt_genrand_res53(&sfmt) - 1); // [0, PI]
                    lig_ori = angle_trans(ph, tht, _radius);
                    dis = dist(lig_ori, ligands.at(j - 1).position_origin);
                } while (dis < clus_minLigDis || dis > clus_maxLigDis);
            }
            else {
                ph = 2.0 * PI * sfmt_genrand_res53(&sfmt);
                tht = acos(2 * sfmt_genrand_res53(&sfmt) - 1); //(0, pi)
                lig_ori = angle_trans(ph, tht, _radius);
            }
        }

        else {
            ph = 2.0 * PI * sfmt_genrand_res53(&sfmt); // [0, 2PI]
            tht = acos(2 * sfmt_genrand_res53(&sfmt) - 1); // [0, PI]
            lig_ori = angle_trans(ph, tht, _radius);
        }

        ligand ligand1;
        ligand1.updatePO(lig_ori, np.position);
        ligands.push_back(ligand1);
    }

}


/**
 * This function binds the NP to the substrate by one bond.
 *
 */
void ini_binding() {

    // Place 1st receptor right underneath the nanoparticle
    receptors.at(0).position.x = np.position.x;
    receptors.at(0).position.y = np.position.y;

    // Place 1st ligand right above the 1st receptor
    ligands.at(0).updatePO(coord{0,0,-_radius}, np.position);

    activeBonds.insert(formBond(0, 0, ligands, receptors, bonds));

    FILE *outfile;
    if ((outfile = fopen(FO7, "w")) == NULL){ printf("\nerror on open FO7!"); exit(0); }
    for (auto j = 0; j < receptorNum; j++)
        fprintf(outfile, "%lf\t%lf\n", receptors.at(j).position.x, receptors.at(j).position.y);
    fclose(outfile);

}

/**
 * This function resumes a simulation from intermediate files
 * "resume.txt", "recepto.txt", "bond_info.txt"
 *
 *
 * @update: timeAcc, sfmt, receptors, ligands, np, bonds, activeBonds
 * @return: if all intermediate files are read correctly
 *
 *
 */
bool resume() {
    double x, y, z;
    int activeBondN, totalBondN;
// read simulation time, particle info, ligands info, active/total bond number
    ifstream input{"resume.txt"};
    if (!input.is_open()) {
        cout << "cant't open resume.txt" << endl;
        return false;
    }
    input >> timeAcc >> np.position.x >> np.position.y >> np.position.z >>
            np.velocity.x >> np.velocity.y >> np.velocity.z >>
            np.rot_velocity.x >> np.rot_velocity.y >> np.rot_velocity.z >> activeBondN >> totalBondN;
    for ( ; !input.eof(); ) {
        input >> x >> y >> z;
        ligand ligand1;
        ligand1.updatePO(coord{x,y,z}, np.position);
        ligands.push_back(ligand1);
    }
    input.close();

// read receptor info
    input.open("receptor.txt");
    if (!input.is_open()) {
        cout << "cant't open receptor.txt" << endl;
        return false;
    }
    for ( ; !input.eof(); ) {
        input >> x >> y;
        receptor receptor1;
        receptor1.position = {x, y, 0};
        receptors.push_back(receptor1);
    }
    input.close();

// read bond info
    input.open("bond_info.txt");
    if (!input.is_open()) {
        cout << "cant't open bond_info.txt" << endl;
        return false;
    }
    int bondNum = 0;
    for( ; !input.eof(); ++bondNum) {
        bond bond;
        input >> bond.bound >> bond.ligand >> bond.receptor >>
                    bond.formPositionLigand.x >> bond.formPositionLigand.y >> bond.formPositionLigand.z >>
                    bond.formPositionReceptor.x >> bond.formPositionReceptor.y >> bond.formPositionReceptor.z >>
                    bond.formTime >> bond.breakTime >>
                    bond.breakPositionLigand.x >> bond.breakPositionLigand.y >> bond.breakPositionLigand.z >>
                    bond.breakPositionReceptor.x >> bond.breakPositionReceptor.y >> bond.breakPositionReceptor.z >>
                    bond.delta;

        bond.name = bondNum;
        bonds.push_back(bond);
        if (bond.bound) {
            receptors.at(bond.receptor).pairing(bond.ligand);
            ligands.at(bond.ligand).pairing(bond.receptor);
            activeBonds.insert(bondNum);
        }
    }
    input.close();

    if ( bonds.size() != totalBondN || activeBonds.size() != activeBondN ) {
        cout << "bond info reads wrong" << endl;
        return false;
    }
    
    if ( receptors.size() != receptorNum || ligands.size() != ligandNum ) {
        cout << "ligand/receptor files read wrong" << endl;
        return false;
    }
    
// assign seed to RNG
    setRNG(sfmt);

    return true;
}

/**
 * This function starts a simulation from one bond state
 *
 * @update: timeAcc, sfmt, receptors, ligands, np, bonds, activeBonds
 *
 */
void ini() {

    timeAcc = 0;
    setRNG(sfmt);
    if (REC_CLU&&DOUBLE_CLU)
        ini_receptor_double_cluster(receptors); // as double clustering
    else if (REC_CLU)
        ini_receptor_single_cluster(receptors); // as single clustering
    else ini_receptor_monomer(receptors); // as monomer

    ini_np(np);
    ini_ligand(ligands);
    ini_binding();
}