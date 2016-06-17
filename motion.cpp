//
// Created by Mingqiu Wang on 5/14/16.
//

#include "declaration.h"

using namespace std;

/**
 * This function calculates the velocity and updates the location of nanoparticle
 *
 * @param: velocity, position, acc
 * @update: velocity, position
 *
 */
void translation(coord & velocity, coord & position, const coord & acc) {

    static const double dev1 = exp(-2.0*beta*timeInc);
    static const double dev2 = (3.0 - 4.0*exp(-beta*timeInc));
    static const double dev3 = 2 - (dev1 + dev2) / (beta*timeInc);
    static const double stddev_pos = sqrt(timeInc*timeInc*thermal / (mass*1e-9) / beta / timeInc*dev3);
    static const double stddev_vel = sqrt(thermal / (mass*1e-9)*(1.0 - exp(-2.0*beta*timeInc))); // (nm/s)
    static const double cor = thermal / (mass*1e-9) / beta * (1 - exp(-beta*timeInc))*
            (1 - exp(-beta*timeInc)) / stddev_pos / stddev_vel; // (dimensionless)
    static const double v_coeff1 = stddev_vel * cor;
    static const double v_coeff2 = stddev_vel * sqrt(1 - cor * cor);
    static const double c0 = exp(-1.0*beta*timeInc); // (dimensionless)
    static const double c1 = (1.0 - c0) / (beta*timeInc); // (dimensionless)
    static const double c2 = (1.0 - c1) / (beta*timeInc); // (dimensionless)
    static const double a = c1 * timeInc;
    static const double b = c2 * timeInc * timeInc * 1e9;
    static const double c = a * 1e9;

    coord ran_x = {gasdev(10000), gasdev(10000), gasdev(10000)};
    coord ran_y = {gasdev(10000), gasdev(10000), gasdev(10000)};

    position = position + a * velocity +  b * acc + ran_x * stddev_pos; // (nm)
    velocity = c0 * velocity + c * acc + v_coeff1 * ran_x + v_coeff2 * ran_y;
}


/**
 * This function calculates the rotational velocity of the particle and generates a rotation matrix
 *
 * @param: velocity, acc
 * @update: velocity
 * @return: rotationMatrix
 */
vector<coord> rotate(coord & velocity, const coord & acc) {

    static const double dev1 = exp(-2.0*beta_rot*timeInc);
    static const double dev2 = (3.0 - 4.0*exp(-beta_rot*timeInc));
    static const double dev3 = 2 - (dev1 + dev2) / (beta_rot*timeInc);
    static const double stddev_pos = sqrt(timeInc*timeInc*thermal / (rot_inertia*1e-9)
                                    / beta_rot / timeInc*dev3);
    static const double stddev_vel = sqrt(thermal / (rot_inertia*1e-9)*(1.0 -
            exp(-2.0*beta_rot*timeInc))); //s^-1
    static const double cor = thermal / (rot_inertia*1e-9) / beta_rot * (1 - exp(-beta_rot*timeInc))*
            (1 - exp(-beta_rot*timeInc)) / stddev_pos / stddev_vel;//dimensional
    static const double v_coeff1 = stddev_vel * cor;
    static const double v_coeff2 = stddev_vel * sqrt(1 - cor * cor);
    static const double c0 = exp(-1.0*beta_rot*timeInc); // (dimensionless)
    static const double c1 = (1.0 - c0) / (beta_rot*timeInc); // (dimensionless)
    static const double c2 = (1.0 - c1) / (beta_rot*timeInc); // (dimensionless)
    static const double a = c1 * timeInc;
    static const double b = c2 * timeInc * timeInc;

    coord ran_x = {gasdev(10000), gasdev(10000), gasdev(10000)};
    coord ran_y = {gasdev(10000), gasdev(10000), gasdev(10000)};

    coord angle = -1 * (a * velocity + b * acc + ran_x * stddev_pos);
    velocity = c0 * velocity + a * acc + v_coeff1 * ran_x + v_coeff2 * ran_y;


// Generates a 3D elementary rotation matrix
    vector<coord> rotationMatrix(3);
    double sinOmega = sin(angle.x);
    double cosOmega = cos(angle.x);
    double sinPhi = sin(angle.y);
    double cosPhi = cos(angle.y);
    double sinKappa = sin(angle.z);
    double cosKappa = cos(angle.z);
    rotationMatrix.at(0).x = cosPhi * cosKappa;
    rotationMatrix.at(0).y = cosOmega * sinKappa + sinOmega * sinPhi * cosKappa;
    rotationMatrix.at(0).z = sinOmega * sinKappa - cosOmega * sinPhi * cosKappa;
    rotationMatrix.at(1).x = -1.0 * cosPhi * sinKappa;
    rotationMatrix.at(1).y = cosOmega * cosKappa - sinOmega * sinPhi * sinKappa;
    rotationMatrix.at(1).z = sinOmega * cosKappa + cosOmega * sinPhi * sinKappa;
    rotationMatrix.at(2).x = sinPhi;
    rotationMatrix.at(2).y = -1.0 * sinOmega * cosPhi;
    rotationMatrix.at(2).z = cosOmega * cosPhi;

    return rotationMatrix;

}

/**
 * This function rotates all ligands on the nanoparticle given a rotation matrix
 *
 * @param: ligands, rotationMatrix, NPposition
 * @update: ligands
 *
 */
void rotateLig (std::vector<ligand> & ligands, const vector<coord> & rotationMatrix, const coord & NPposition) {
    double mol_cal[3];
    for (auto & ligand: ligands) {
        for (auto ii = 0; ii < 3; ii++) {
            mol_cal[ii] =
                    rotationMatrix.at(ii).x * ligand.position_origin.x +
                    rotationMatrix.at(ii).y * ligand.position_origin.y +
                    rotationMatrix.at(ii).z * ligand.position_origin.z;
        }
        ligand.updatePO({mol_cal[0], mol_cal[1], mol_cal[2]}, NPposition);
    }
}