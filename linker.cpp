//
// Created by Mingqiu Wang on 6/21/16.
//


#include "declaration.h"

/**
 * This
 *
 * @update:
 *
 */
double delta2_com(double delta) {
    if (_maxDeltaCom + delta > 0 ) return 0;
    return delta + _maxDeltaCom;
}


double delta2_ext(double delta) {
    double delta1 = _sigma * delta / (_sigma + _sigma_linker);
    if (delta1 > _maxDeltaExt) return delta - _maxDeltaExt; // assume the linker length is far below bond length
    return _sigma_linker * delta / (_sigma + _sigma_linker);
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
        if (ifBreak(deltaBond)) {
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

/**
 * This function calculates force from all bonds
 *
 * @param:
 *
 * @return: pair<coord, coord>, 1st is force, 2nd is torque
 *
 */
std::pair<coord, coord> Fbond_linker(const std::set<int> & activeBond, const std::vector<bond> & bonds,
                              const std::vector<receptor> & receptors, const std::vector<ligand> & ligands) {

    int lig, rec;
    double rot_f_x, rot_f_y, rot_f_z;
    coord force, sumF, sumT;

    for (int bond : activeBond) {
        lig = bonds.at(bond).ligand;
        rec = bonds.at(bond).receptor;
        force = _sigma * bonds.at(bond).delta / dist(receptors.at(rec).position, ligands.at(lig).position) *
                (receptors.at(rec).position - ligands.at(lig).position); // (nN)
        sumF = sumF + force;
        rot_f_x = (ligands.at(lig).position_origin.y * force.z - ligands.at(lig).position_origin.z * force.y); // (nN*nm)
        rot_f_y = (ligands.at(lig).position_origin.z * force.x - ligands.at(lig).position_origin.x * force.z); // (nN*nm)
        rot_f_z = (ligands.at(lig).position_origin.x * force.y - ligands.at(lig).position_origin.y * force.x); // (nN*nm)
        // right-hand system
        sumT = sumT + coord{rot_f_x, rot_f_y, rot_f_z}; // (nN*nm)
    }
    return {sumF, sumT};
}