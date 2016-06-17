//
// Created by Mingqiu Wang on 5/7/16.
//


#include "declaration.h"

/**
 * This function sets up all necessities for the simulation:
 *      1) prepares the seed for random number generator: sfmt
 *      2) defines the bonding length truncation: bondCutoff
 *      3) sets system time: timeAcc
 *
 * If in DEBUG mode,
 *      1) overwrites output files
 *      2) sets seed for random number generator as 1
 *
 * @update: timeAcc, sfmt, bondCutoff
 * @file: seed.txt
 *
 */
void setup() {

    FILE *outfile;
    timeAcc = 0;

    if (DEBUG) {
        if ((outfile = fopen(FO1, "w")) == NULL){ printf("\nerror on open FO1!"); exit(1); }
        fclose(outfile);
        sfmt_init_gen_rand(&sfmt, 1);
    }

    else {
        unsigned long long seed;
	    FILE *fp = fopen("/dev/urandom", "r"); fread(&seed, 1, sizeof(unsigned int), fp); fclose(fp);
	    if ((outfile = fopen(FO3, "w")) == NULL){ printf("\nerror on open FO3!"); exit(1); }
	    fprintf(outfile, "%e\t%llu\n", timeAcc, seed);
	    fclose(outfile);
        sfmt_init_gen_rand(&sfmt, seed);
    }

    bondCutoff.deltaMax = 0.05 * bondL;
    bondCutoff.bondLMax = bondL + bondCutoff.deltaMax;
    bondCutoff.bondLMin = bondL - bondCutoff.deltaMax; //TODO: a better cutoff

}


/**
 * This function distributes a point in 2D close to a given point i
 * (closer than upper_dis && further than lower_dis).
 *
 * @param
 * @param
 * @return
 */
std::pair<double, double> distribute
        (double upper_dis, double lower_dis, const std::pair<double, double>& i) {
    double x, y, d;
    do {
        x = sfmt_genrand_res53(&sfmt) * 2 * upper_dis + i.first - upper_dis;
        y = sfmt_genrand_res53(&sfmt) * 2 * upper_dis + i.second - upper_dis;
        d = dist(x, y, 0, i.first, i.second, 0);

    } while ((d>upper_dis) || (d<lower_dis));
    // (x, y) and i is further away from upper_dis or closer than lower_dis

    return{ x, y };

}


/**
 * This function calculates the distance between two points.
 *
 * @param
 * @param
 * @return
 */
double dist(double rectx, double recty, double rectz, double surfx, double surfy, double surfz) {
    return sqrt((rectx - surfx)*(rectx - surfx) + (recty - surfy)*(recty - surfy) + (rectz - surfz)*(rectz - surfz));
}


/**
 * This function calculates the distance between two points.
 *
 * @param
 * @param
 * @return
 */
double dist(const coord &coord1, const coord &coord2) {
    return sqrt((coord1.x - coord2.x)*(coord1.x - coord2.x) +
                (coord1.y - coord2.y)*(coord1.y - coord2.y) +
                (coord1.z - coord2.z)*(coord1.z - coord2.z));
}


/**
 * This function transfers spherical coordinate (ph, tht, radius) to
 * Cartesian coordinate (coord)
 * @param: ph [0, 2PI], tht [0, PI], radius
 * @return
 */
coord angle_trans(double ph, double tht, double radius) {
    return { radius*sin(tht)*cos(ph), radius*sin(tht)*sin(ph), radius*cos(tht)};
}

/**
 * This function collects all available receptors within
 * radius + 2 * bond length of NP
 *
 * @update: availRec, np.lastPairPos
 */
void getAvailRec() {
    availRec.clear();
    for (auto j = 0; j < receptorNum; ++j)
        if (dist(np.position,receptors.at(j).position)<(radius + 2 * bondL))
            availRec.push_back(j);
    np.lastPairPos = np.position;
}

/**
 * This function collects all available ligands below bond length
 *
 * @update: availLig
 */
void getAvailLig() {
    availLig.clear();
    for (auto j = 0; j < ligandNum; ++j)
        if (ligands.at(j).position.z<bondL)
            availLig.push_back(j);
}


/**
 * File:      gasdev.C
 * Date:      12/31/01  (written fall 2001)
 * Author:    C. Paciorek
 * Purpose:   generate a standard normally-distributed (i.e. Gaussian) random variable (scalar)
 * Reference: Numerical Recipes in C - code taken directly from the source
 */
double gasdev(long idum) {  // generates 2 random deviates at a time; idum is seed
    static int iset=0;
    static double gset;
    double fac,rsq,v1,v2;

    if (idum < 0) iset=0;
    if (iset == 0) {  //We don't have an extra deviate handy, so
        do {
            v1=2.0*sfmt_genrand_res53(&sfmt)-1.0; // pick two uniform numbers in the square extending from -1 to +1 in each direction,
            v2=2.0*sfmt_genrand_res53(&sfmt)-1.0;
            rsq=v1*v1+v2*v2; // see if they are in the unit circle,
        }
        while (rsq >= 1.0 || rsq == 0.0); // and if they are not, try again.
        fac=sqrt(-2.0*log(rsq)/rsq); /* Now make the Box-Muller transformation to get two normal deviates. Return one and save the other for next time. */
        gset=v1*fac;
        iset=1; // Set ag.
        return v2*fac;
    }
    else { // We have an extra deviate handy,
        iset=0; // so unset the ag,
        return gset; // and return it.
    }
}