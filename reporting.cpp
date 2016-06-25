//
// Created by Mingqiu Wang on 5/26/16.
//

#include "declaration.h"

using namespace std;
bool checkDisplace(long &step) {

    if (dist(np.lastPairPos,np.position) > _bondEL) getAvailRec(); // get all available receptors
    getAvailLig(); // get all available ligands

    if (ifDetach(np.position)) return true;

    if (!(step%TRAJ)) {
        writeLoc();
    }
    if (!(step%REPORTER)) {
        writeBond();
        writeResume();
    }
    return false;

}

void writeLoc() {

    FILE *outfile;
    if ((outfile = fopen(FO1, "a")) == NULL){ printf("\nerror on open FO1!"); exit(0); }
    fprintf(outfile, "%.4e\t%lf\t%lf\t%lf\t%lu\n",
            timeAcc, np.position.x, np.position.y, np.position.z, activeBond.size());
    fclose(outfile);

}


void writeEndTime() {

    FILE *outfile;
    if ((outfile = fopen(FO5, "w")) == NULL){ printf("\nerror on open FO5!"); exit(0); }
    fprintf(outfile, "%.4e\t%lf\t%lf\t%lf\t%lu\n", timeAcc, np.position.x,
            np.position.y, np.position.z, activeBond.size());
    fclose(outfile);

}


/**
 * This function writes bond information into the file
 *
 * @file: bond_info.txt
 *      (0th column: bound status
 *      1st column: ligand
 *      2nd column: receptor
 *      3rd/4th/5th columns: ligand coordinate when bond formed
 *      6th/7th/8th columns: receptor coordinate when bond formed
 *      9th column: time when bond formed
 *      10th column: time when bond broke (-1 if the bond hasn't broken until the end of simulation)
 *      11st/12nd/13rd columns: ligand coordinate when bond broke
 *      14th/15th/16th columns: receptor coordinate when bond broke)
 *
 */
void writeBond() {

    FILE *outfile;
    if ((outfile = fopen(FO6, "w")) == NULL){ printf("\nerror on open FO6!"); exit(0); }
    for (const auto & bond : bonds)
        fprintf(outfile, "%d\t%d\t%d\t"
                        "%lf\t%lf\t%lf\t"
                        "%lf\t%lf\t%lf\t"
                        "%.4e\t%.4e\t"
                        "%lf\t%lf\t%lf\t"
                        "%lf\t%lf\t%lf\t"
                        "%lf\n",
                bond.bound, bond.ligand, bond.receptor,
                bond.formPositionLigand.x, bond.formPositionLigand.y, bond.formPositionLigand.z,
                bond.formPositionReceptor.x, bond.formPositionReceptor.y, bond.formPositionReceptor.z,
                bond.formTime, bond.breakTime,
                bond.breakPositionLigand.x, bond.breakPositionLigand.y, bond.breakPositionLigand.z,
                bond.breakPositionReceptor.x, bond.breakPositionReceptor.y, bond.breakPositionReceptor.z,
                bond.delta
        );

    fclose(outfile);

}

void writeResume() {

    FILE *outfile;
    if ((outfile = fopen(FO8, "w")) == NULL){ printf("\nerror on open FO8!"); exit(0); }

    fprintf(outfile, "%.4e\n%lf\t%lf\t%lf\n"
                        "%lf\t%lf\t%lf\n"
                        "%lf\t%lf\t%lf\n",
                timeAcc, np.position.x, np.position.y, np.position.z,
                np.velocity.x, np.velocity.y, np.velocity.z,
                np.rot_velocity.x, np.rot_velocity.y, np.rot_velocity.z);
    for (const auto & ligand : ligands)
        fprintf(outfile, "%lf\t%lf\t%lf\n",
                ligand.position_origin.x, ligand.position_origin.y,
                ligand.position_origin.z);


    fclose(outfile);

}