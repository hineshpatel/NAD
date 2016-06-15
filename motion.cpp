//
// Created by Mingqiu Wang on 5/14/16.
//

#include "declaration.h"

using namespace std;

/**
 * This function updates the location of nanoparticle
 *
 * @coefficient: c0, c1, c2, v_coeff1, v_coeff2, stddev_pos
 * @param[global var]: np.position (old), np.velocity (old)
 * @return[global var]: np.position (new), np.velocity (new)
 *
 */
void translate() {

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

    np.position = np.position + a * np.velocity +  b * np.acc + ran_x * stddev_pos; // (nm)
    np.velocity = c0 * np.velocity + c * np.acc + v_coeff1 * ran_x + v_coeff2 * ran_y;
}


/**
 * This function rotates the nanoparticle, and also updates the location of each ligand
 *
 * @coefficient: c0, c1, c2, v_coeff1, v_coeff2, stddev_pos
 * @param[global var]: np.angle (old), np.rot_velocity (old), ligands(old)
 * @return[global var]: np.angle (new), np.rot_velocity (new), ligands(old)
 *
 */
void rotate() {

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

    np.angle = -1 * (a * np.rot_velocity + b * np.rot_acc + ran_x * stddev_pos);
    np.rot_velocity = c0 * np.rot_velocity + a * np.rot_acc + v_coeff1 * ran_x + v_coeff2 * ran_y;


// updates 3D elementary rotation matrix
    vector<coord> rotationMatrix(3);
    double mol_cal[3];
    double sinOmega = sin(np.angle.x);
    double cosOmega = cos(np.angle.x);
    double sinPhi = sin(np.angle.y);
    double cosPhi = cos(np.angle.y);
    double sinKappa = sin(np.angle.z);
    double cosKappa = cos(np.angle.z);
    rotationMatrix.at(0).x = cosPhi * cosKappa;
    rotationMatrix.at(0).y = cosOmega * sinKappa + sinOmega * sinPhi * cosKappa;
    rotationMatrix.at(0).z = sinOmega * sinKappa - cosOmega * sinPhi * cosKappa;
    rotationMatrix.at(1).x = -1.0 * cosPhi * sinKappa;
    rotationMatrix.at(1).y = cosOmega * cosKappa - sinOmega * sinPhi * sinKappa;
    rotationMatrix.at(1).z = sinOmega * cosKappa + cosOmega * sinPhi * sinKappa;
    rotationMatrix.at(2).x = sinPhi;
    rotationMatrix.at(2).y = -1.0 * sinOmega * cosPhi;
    rotationMatrix.at(2).z = cosOmega * cosPhi;

    for (auto & ligand: ligands) {
        for (auto ii = 0; ii < 3; ii++) {
            mol_cal[ii] =
            rotationMatrix.at(ii).x * ligand.position_origin.x +
            rotationMatrix.at(ii).y * ligand.position_origin.y +
            rotationMatrix.at(ii).z * ligand.position_origin.z;
        }
        ligand.updatePO({mol_cal[0], mol_cal[1], mol_cal[2]}, np.position);
    }

}