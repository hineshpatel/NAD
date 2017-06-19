//
// Created by Mingqiu Wang on 5/14/16.
//

#include "declaration.h"


/**
 * This function checks all active bonds and determines possible breakage.
 *
 * @update: activeBonds; ligands; receptors; bonds
 *
 */
void breakageCheck(std::set<int> & activeBonds, std::vector<bond> & bonds,
                   std::vector<ligand> & ligands, std::vector<receptor> & receptors) {
    double delta;
    int lig, rec;

    for (auto bond = activeBonds.begin(); bond != activeBonds.end();) {
        lig = bonds.at(*bond).ligand;
        rec = bonds.at(*bond).receptor;
        coord ligStem = (ligands.at(lig).position - np.position) * (_radius /
                dist(ligands.at(lig).position, np.position)) + np.position;
        delta = dist(ligStem, receptors.at(rec).stem) - _bondEL; // (nm)
        bonds.at(*bond).delta = delta;
        if (ORI) {
            if (delta < 0) {
                ++bond;
                continue;
            }
        }
        if (ifBreak(delta)) {
            ligands.at(lig).unpairing();
            receptors.at(rec).unpairing();
            bonds.at(*bond).bound = false;
            bonds.at(*bond).breakTime = timeAcc;
            bonds.at(*bond).breakPositionLigand = ligStem;
            bonds.at(*bond).breakPositionReceptor = receptors.at(rec).stem;
            bonds.at(*bond).delta = -1.0;
            bond = activeBonds.erase(bond);
        }
        else ++bond;
    }
}



/**
 * This function determines if the bond breaks given the deviation from equilibrium length
 * Equation: kr = kr0 * exp (gamma * sigma * delta / kB / T)
 *
 * @param: delta: the length that the bond is stretched or compressed
 *      from the equilibrium length
 * @return: if breaks
 *
 */
bool ifBreak(double delta) {
    double kr;
    static const double kr_cal = (_sigma * 1000.0 * _gama) / _thermal; // (nm^-1)

    kr = _kr0 * exp(kr_cal * fabs(delta));
    return (sfmt_genrand_res53(&sfmt)< (1.0 - exp(-1.0 * kr * _timeInc)));

}



/**
 * This function checks available unbound adhesion molecules and
 *      determines possible formation.
 *
 * @param: availLig: available ligands; availRec: available receptors
 * @update: activeBonds; ligands; receptors
 *
 */
void formationCheck(const std::vector<int> & availLig, const std::vector<int> & availRec,
                    std::set<int> & activeBonds, std::vector<ligand> & ligands,
                    std::vector<receptor> & receptors) {

    double bondDis;

    for (auto lig: availLig) {
        if (ligands.at(lig).bound) continue;
        for (auto rec: availRec) {
            if (receptors.at(rec).bound) continue;
            bondDis = dist(ligands.at(lig).position,
                           receptors.at(rec).position); // (nm)
            if (bondDis > bondCutoff.bondLMax || bondDis < bondCutoff.bondLMin) continue;
            if (ifForm(bondDis-_bondEL)) {
                if (CROSSCHECK&&ifCross (activeBonds, receptors, ligands, lig, rec)) continue;
                activeBonds.insert(formBond(lig, rec, ligands, receptors, bonds));
                break;
            }
        }
    }
}

/**
 * This function checks available unbound adhesion molecules and
 *      determines possible formation for Ori enabled condition only.
 *
 * @param: availLig: available ligands; availRec: available receptors
 * @update: activeBonds; ligands; receptors
 *
 */
void formationCheckOri(const std::vector<int> & availLig, const std::vector<int> & availRec,
                    std::set<int> & activeBonds, std::vector<ligand> & ligands,
                    std::vector<receptor> & receptors) {

    double bondDis;

    for (auto lig: availLig) {
        if (ligands.at(lig).bound) continue;
        for (auto rec: availRec) {
            if (receptors.at(rec).bound) continue;
            bondDis = dist(ligands.at(lig).position,
                           receptors.at(rec).position); // (nm)
            if (bondDis > bondCutoff.deltaMax || bondDis < -1.0*bondCutoff.deltaMax) continue;
            if (ifForm(bondDis)) {
                if (CROSSCHECK&&ifCrossOri (activeBonds, receptors, ligands, lig, rec)) continue;
                activeBonds.insert(formBond(lig, rec, ligands, receptors, bonds));
                break;
            }
        }
    }
}

/**
 * This function determines if a potential bond would form
 *      given the distance between pairing ligand and receptor
 * Equation: kf = kf0 * exp ( - sigma_ts * delta * delta / (2*kB*T))
 *
 * @param: bondLength: the distance between pairing ligand and receptor
 * @return: if forms
 *
 */
bool ifForm(double delta) {
    double kf;
    static const double kf_cal = -1.0 * (sigma_ts * 500.0) / _thermal; // (nm^-2)

    kf = _kf0 * exp(kf_cal * delta * delta);
    if (sfmt_genrand_res53(&sfmt)< 1.0 - exp(-1.0 * kf * _timeInc)) {
        return !ifBreak(delta);
    }
    return false;
}

/**
 * This function forms a bond between receptor rec and ligand lig.
 *
 * @param: lig: No. lig and rec: No. rec
 * @return: No. the bond
 * @update: vector<ligand> ligands, vector<receptor> receptors, vector<bond> bonds
 *
 */
int formBond(int lig, int rec, std::vector<ligand> & ligands, std::vector<receptor> & receptors,
             std::vector<bond> & bonds) {

    ligands.at(lig).pairing(rec);
    receptors.at(rec).pairing(lig);

    coord ligStem = (ligands.at(lig).position - np.position) * (_radius /
          dist(ligands.at(lig).position, np.position)) + np.position;

    bond bond1;
    bond1.name = bonds.size();
    bond1.bound = true;
    bond1.formPositionLigand = ligStem;
    bond1.formPositionReceptor = receptors.at(rec).stem;
    bond1.formTime = timeAcc;
    bond1.ligand = lig;
    bond1.receptor = rec;
    bond1.delta = dist(bond1.formPositionLigand, bond1.formPositionReceptor) - _bondEL;
    bonds.push_back(bond1);

    return bond1.name;

}


/**
 * This function checks if the particle moves outside the substrate
 * assume the initial position of NP is at (x = 0, y = 0)
 *
 * @param: position: current position of NP
 * @return: if the particle moves more than (_detach_criteria) radius
 */
bool ifDetach(const coord & position) {
    const static double dis = _detach_criteria * _radius;
    return (dist({position.x, position.y}, {0, 0}) > dis);
}
