//
// Created by Mingqiu Wang on 5/7/16.
//

#include "declaration.h"

/**
 * This function initializes receptors on the substrate under double clustering condition.
 *
 * Initializes: std::vector<receptor> receptors
 * Outputs: receptor.txt
 *
 */
void ini_receptor_double_cluster() {

    int gProteinNum = gProteinDens * substrate * substrate / 250000; // number of G protein sites
    std::vector<int> gProtein (gProteinNum + 1, 0); // how many receptors in each G protein site (max 4)
    std::vector<std::pair<double, double> > gProteinPos; // positions of G protein sites (nm)
    double x, y; srand(time(NULL));

    // distributes G protein sites on the substrate
    for (auto j = 0; j < gProteinNum; j++) {
        x = sfmt_genrand_res53(&sfmt) * 2* substrate - substrate;
        y = sfmt_genrand_res53(&sfmt) * 2* substrate - substrate;
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
 * Initializes: std::vector<receptor> receptors
 * Outputs: receptor.txt
 *
 */
void ini_receptor_single_cluster() {

    double receptor_x, receptor_y;
    for (auto j = 0; j < receptorNum; j++) { // distribute receptors
        if (j % 2) {
            std::pair<double, double> a = distribute(clus_maxRecDis,
                                                     clus_minRecDis, {receptors.at(j - 1).position.x, receptors.at(j - 1).position.y });
            receptor_x = a.first;
            receptor_y = a.second;
        }
        else {
            receptor_x = sfmt_genrand_res53(&sfmt) * 2* substrate - substrate;
            receptor_y = sfmt_genrand_res53(&sfmt) * 2* substrate - substrate;
        }
        receptor receptor1;
        receptor1.position = {receptor_x, receptor_y, 0};
        receptors.push_back(receptor1);
    }
}


/**
 * This function initializes receptors on the substrate under monomer condition.
 *
 * Initializes: std::vector<receptor> receptors
 * Outputs: receptor.txt
 *
 */
void ini_receptor_monomer() {

    double receptor_x, receptor_y;
    for (auto j = 0; j < receptorNum; j++) { // distribute receptors

        receptor_x = sfmt_genrand_res53(&sfmt) * 2* substrate - substrate;
        receptor_y = sfmt_genrand_res53(&sfmt) * 2* substrate - substrate;

        receptor receptor1;
        receptor1.position = {receptor_x, receptor_y, 0};
        receptors.push_back(receptor1);
    }

}


/**
 * This function initializes a nanoparticle.
 *
 * Initializes: np np
 *
 */
void ini_np() {

    np.position.x = 0;
    np.position.y = 0;
    np.position.z = radius + bondL;
    np.lastPairPos = np.position; //TODO: is this last pair pos useful?

}


/**
 * This function initializes ligandss on the NP.
 *
 * Initializes: std::vector<ligand> ligands
 *
 */
void ini_ligand() {

    double ph, tht, dis; coord lig_ori;
    for (auto j = 0; j < ligandNum; j++) {
        if (LIG_CLU) {
            if (j % 2) {
                do {
                    ph = 2.0 * PI * sfmt_genrand_res53(&sfmt); // [0, 2PI]
                    tht = acos(2 * sfmt_genrand_res53(&sfmt) - 1); // [0, PI]
                    lig_ori = angle_trans(ph, tht, radius);
                    dis = dist(lig_ori, ligands.at(j - 1).position_origin);
                } while (dis < clus_minLigDis || dis > clus_maxLigDis);
            }
            else {
                ph = 2.0 * PI * sfmt_genrand_res53(&sfmt);
                tht = acos(2 * sfmt_genrand_res53(&sfmt) - 1); //(0, pi)
                lig_ori = angle_trans(ph, tht, radius);
            }
        }

        else {
            ph = 2.0 * PI * sfmt_genrand_res53(&sfmt); // [0, 2PI]
            tht = acos(2 * sfmt_genrand_res53(&sfmt) - 1); // [0, PI]
            lig_ori = angle_trans(ph, tht, radius);
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
    ligands.at(0).updatePO(coord{0,0,-radius}, np.position);

    activeBond.insert(formBond(0, 0));

    FILE *outfile;
    if ((outfile = fopen(FO7, "w")) == NULL){ printf("\nerror on open FO7!"); exit(0); }
    for (auto j = 0; j < receptorNum; j++)
        fprintf(outfile, "%lf\t%lf\n", receptors.at(j).position.x, receptors.at(j).position.y);
    fclose(outfile);

}
