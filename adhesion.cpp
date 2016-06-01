//
// Created by Mingqiu Wang on 5/14/16.
//

#include "declaration.h"

void breakageCheck() {
    double bondDis, kr, prob;
    int lig, rec;
    static double kr_cal = (sigma * 1000.0 * gama) / thermal; // (nm^-2)

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
            bonds.at(*bond).bind = false;
            bonds.at(*bond).breakTime = timeAcc;
            bonds.at(*bond).breakPositionLigand = ligands.at(lig).position;
            bonds.at(*bond).breakPositionReceptor = receptors.at(rec).position;
            bond = activeBond.erase(bond);

        }
        else ++bond;
    }

}

void formationCheck() {

    double bond_distance, kf, kr, prob;
    static double kr_cal = (sigma * 1000.0 * gama) / thermal; // (nm^-2)
    const double kf_cal = -1.0*(sigma_ts*500.0) / thermal; // (nm^-2)

    for (auto ligand: availLig) {
        if (ligands.at(ligand).bind) continue;
        for (auto receptor: availRec) {
            if (receptors.at(receptor).bind) continue;
            bond_distance = dist(ligands.at(ligand).position,
            receptors.at(receptor).position); // (nm)
            kf = kf0*exp(kf_cal*((bond_distance - bondL)*(bond_distance - bondL)));
            prob = 1.0 - exp(-1.0*kf*timeInc);
            if (sfmt_genrand_res53(&sfmt)< prob) {
                kr = kr0*exp(kr_cal*fabs(bond_distance - bondL));
                prob = 1.0 - exp(-1.0*kr*timeInc);
                if (sfmt_genrand_res53(&sfmt)> prob) {
                    activeBond.insert(formBond(ligand, receptor));
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
 *
 */
int formBond(int lig, int rec) {

    ligands.at(lig).pairing(rec);
    receptors.at(rec).pairing(lig);

    bond bond1;
    bond1.name = bonds.size();
    bond1.bind = true;
    bond1.formPositionLigand = ligands.at(lig).position;
    bond1.formPositionReceptor = receptors.at(rec).position;
    bond1.formTime = timeAcc;
    bond1.ligand = lig;
    bond1.receptor = rec;
    bonds.push_back(bond1);

    return bond1.name;

}

