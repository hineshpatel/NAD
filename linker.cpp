//
// Created by Mingqiu Wang on 6/21/16.
//


#include "declaration.h"

/**
 * This function checks all active bonds and determines possible breakage.
 * Equation: kr = kr0 * exp (gamma * sigma * delta / kB / T)
 *
 * @update:
 *
 */
double delta2_com(double delta) {
    delta = fabs(delta);
    if (delta < maxDeltaCom) return 0;
    return - delta + maxDeltaCom;
}


double delta2_ext(double delta) {
    delta = fabs(delta);
    double delta1 = _sigma * delta / (_sigma + sigma_linker);
    if (delta1 > maxDeltaExt) return delta - maxDeltaExt; // assume the linker length is far below bond length
    return sigma_linker * delta / (_sigma + sigma_linker);
}

/**
 * This function checks all active bonds and determines possible breakage.
 * Equation: kr = kr0 * exp (gamma * sigma * delta / kB / T)
 *
 * @update:
 *
 */
void breakageCheck_linker(std::set<int> & activeBonds, std::vector<bond> & bonds,
                          std::vector<ligand> & ligands, std::vector<receptor> & receptors) {

    double delta, deltaBond;
    int lig, rec;

    for (auto bond = activeBonds.begin(); bond != activeBonds.end();) {
        lig = bonds.at(*bond).ligand;
        rec = bonds.at(*bond).receptor;
        delta = dist(ligands.at(lig).position, receptors.at(rec).position) - _bondEL; // (nm)
        if (delta < 0) deltaBond = delta2_com(delta);
        else deltaBond = delta2_ext(delta);
        bonds.at(*bond).delta = deltaBond;
        if (ifBreak(fabs(deltaBond))) {
            ligands.at(lig).unpairing();
            receptors.at(rec).unpairing();
            bonds.at(*bond).bound = false;
            bonds.at(*bond).breakTime = timeAcc;
            bonds.at(*bond).breakPositionLigand = ligands.at(lig).position;
            bonds.at(*bond).breakPositionReceptor = receptors.at(rec).position;
            bonds.at(*bond).delta = -1.0;
            bond = activeBonds.erase(bond);
        }
        else ++bond;
    }
}