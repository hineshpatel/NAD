//
// Created by Mingqiu Wang on 5/14/16.
//

#include "declaration.h"


/**
 * This function checks all active bonds and determines possible breakage.
 * kr = kr0 * exp (gamma * sigma * delta / kB / T)
 *
 * @coefficient: kr_cal
 * @update: set<int> activeBond, vector<ligand> ligands,
 *      vector<receptor> receptors, vector<bond> bonds
 *
 */
void breakageCheck() {
    double bondDis, kr, prob;
    int lig, rec;
    static const double kr_cal = (sigma * 1000.0 * gama) / thermal; // (nm^-1)

    for (auto bond = activeBond.begin(); bond != activeBond.end();) {
        lig = bonds.at(*bond).ligand;
        rec = bonds.at(*bond).receptor;
        bondDis = dist(ligands.at(lig).position,
                       receptors.at(rec).position); // (nm)
        kr = kr0*exp(kr_cal*fabs(bondDis - bondL));
        prob = 1.0 - exp(-1.0*kr*timeInc);
        if (sfmt_genrand_res53(&sfmt)< prob) {
            ligands.at(lig).unpairing();
            receptors.at(rec).unpairing();
            bonds.at(*bond).bound = false;
            bonds.at(*bond).breakTime = timeAcc;
            bonds.at(*bond).breakPositionLigand = ligands.at(lig).position;
            bonds.at(*bond).breakPositionReceptor = receptors.at(rec).position;
            bond = activeBond.erase(bond);
        }
        else ++bond;
    }
}

/**
 * This function checks available unbound adhesion molecules and
 *      determines possible formation.
 * kf = kf0 * exp ( - sigma_ts * delta * delta / (2*kB*T))
 *
 * @coefficient: kf_cal, kr_cal
 * @update: set<int> activeBond, vector<ligand> ligands,
 *      vector<receptor> receptors, vector<bond> bonds
 *
 */
void formationCheck() {

    double bondDis, kf, kr, prob;
    static const double kr_cal = (sigma * 1000.0 * gama) / thermal; // (nm^-1)
    static const double kf_cal = -1.0 * (sigma_ts * 500.0) / thermal; // (nm^-2)

    for (auto lig: availLig) {
        if (ligands.at(lig).bound) continue;
        for (auto rec: availRec) {
            if (receptors.at(rec).bound) continue;
            bondDis = dist(ligands.at(lig).position,
                           receptors.at(rec).position); // (nm)
            kf = kf0*exp(kf_cal*((bondDis - bondL) * (bondDis - bondL)));
            prob = 1.0 - exp(-1.0*kf*timeInc);
            if (sfmt_genrand_res53(&sfmt)< prob) {
                kr = kr0*exp(kr_cal*fabs(bondDis - bondL));
                prob = 1.0 - exp(-1.0*kr*timeInc);
                if (sfmt_genrand_res53(&sfmt)> prob) {
                    activeBond.insert(formBond(lig, rec));
                    break;
                }
            }
        }
    }
}


/**
 * This function forms a bond between receptor rec and ligand lig.
 *
 * @param: no. lig and no. rec
 * @return: no. the bond
 * @update: vector<ligand> ligands, vector<receptor> receptors, vector<bond> bonds
 *
 */
int formBond(int lig, int rec) {

    ligands.at(lig).pairing(rec);
    receptors.at(rec).pairing(lig);

    bond bond1;
    bond1.name = bonds.size();
    bond1.bound = true;
    bond1.formPositionLigand = ligands.at(lig).position;
    bond1.formPositionReceptor = receptors.at(rec).position;
    bond1.formTime = timeAcc;
    bond1.ligand = lig;
    bond1.receptor = rec;
    bonds.push_back(bond1);

    return bond1.name;

}


/**
 * This function checks if the particle moves outside the substrate
 *
 */
bool ifDetach() {
    const static double dis = detach_criteria * radius;
    return (dist(np.position, coord{}) > dis);
}
