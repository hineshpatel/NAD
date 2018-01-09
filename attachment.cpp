//
// Created by Mingqiu Wang on 11/2/17.
//


#include <vector>
#include "dataStructures.h"
#include "parameters.h"
#include "declaration.h"

/**
 * This function checks whether NP is still in the initial cell.
 *
 */
bool inCell(const coord& NPposition) {
    return NPposition.x<_boxLength && NPposition.x>-1*_boxLength &&
    NPposition.y<_boxLength && NPposition.y>-1*_boxLength;
}

/**
 * This function puts the NP back to its box.
 */
void putNPBack(coord& NPposition) {
    int adjust = int(NPposition.x/_boxLength);
    NPposition.x -= adjust*_boxLength*2;
    if (!(adjust == 1 || adjust == -1 || adjust == 0))
        putNPBack(NPposition);
    adjust = int(NPposition.y/_boxLength);
    NPposition.y -= adjust*_boxLength*2;
    if (!(adjust == 1 || adjust == -1 || adjust == 0))
        putNPBack(NPposition);
}

/**
 * This function adds a new NP into the box.
 */
void renewNP() {
    // set all receptors unbind, so far no need
    ini_np_rand(np);
    ini_ligand(ligands);
    // get available adhesion molecules
    getAvailRec(availRec, np);
    getAvailLig(availLig, ligands);
}



void writeAttNum(int num) {

    FILE *outfile;
    if ((outfile = fopen(FO11, "a")) == NULL){ printf("\nerror on open FO11!"); exit(0); }
    fprintf(outfile, "%.4e\t%d\t%d\n", timeAcc, num, np.name);
    fclose(outfile);
    printf("%.4e\t%d\t%d\n", timeAcc, num, np.name);

}

void summarizeNP(double start_t, double end_t, bool attached) {

    FILE *outfile;
    if ((outfile = fopen(FO12, "a")) == NULL){ printf("\nerror on open FO12!"); exit(0); }
    fprintf(outfile, "%.4e\t%.4e\t%d\n", start_t, end_t, attached);
    fclose(outfile);
}