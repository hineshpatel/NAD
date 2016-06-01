//
// Created by Mingqiu Wang on 5/14/16.
//

#include "declaration.h"

/**
 * This function calculates the repulsion force
 *
 * @param: np.position.z
 * @return: repulsion force (nN)
 *
 */
coord Frepulsion() {

    double distance = (np.position.z - radius);
    if (distance <= bondL)
    {
        double contact_area = PI*(radius*radius - (np.position.z - bondL)*(np.position.z - bondL));
        double force = compressibility*(1 / distance*(1 / distance + 1/polymer_thickness))*
                exp(-1.0*distance*1/polymer_thickness);
        //r*exp(-s/t)*(1/s/s+1/t/s)
        return coord{0, 0, contact_area*force};
    }
    else return coord{};

}

/**
 * This function calculates the bond force
 *
 * @param:
 *
 * @return: pair<coord, coord>, 1st is force, 2nd is torque
 *
 */
std::pair<coord, coord> Fbond() {

    int lig, rec;
    double xb, force_cal, rot_f_x, rot_f_y, rot_f_z;
    coord force, torque;
    for (int bond : activeBond) {
        lig = bonds.at(bond).ligand;
        rec = bonds.at(bond).receptor;
        xb = dist(ligands.at(lig).position, receptors.at(rec).position); // (nm)
        force = force + sigma*(xb - bondL) / xb *
                (receptors.at(rec).position - ligands.at(lig).position); //(nN)
        rot_f_x =  (ligands.at(lig).position_origin.y * force.z - (ligands.at(lig).position.z) * force.y); // (nN*nm)
        rot_f_y = (ligands.at(lig).position_origin.z * force.x - (ligands.at(lig).position.x) * force.z); // (nN*nm)
        rot_f_z = (ligands.at(lig).position_origin.x * force.y - (ligands.at(lig).position.y) * force.x); // (nN*nm)
        // right-hand system
        torque = torque + coord{rot_f_x, rot_f_y, rot_f_z}; // (nN*nm)
    }
    return {force, torque};

}


std::pair<coord, coord> Fshear() {

    static double trans_shear_coef = 0.006*PI*viscosity*radius*shear_rate; // nN/nm
    static double trans_shear_coef2 = 0.006*0.5625*PI*viscosity*radius*radius*shear_rate; // nN/nm
    static double rot_shear_coef = 0.004 * PI*viscosity*radius*radius*radius*shear_rate; // nN*nm
    double trans_shear;
    double rot_shear, rot_shear_coef2;
    trans_shear = trans_shear_coef* np.position.z + trans_shear_coef2;
    rot_shear_coef2 = radius / np.position.z;
    rot_shear = rot_shear_coef*(1 - 0.1875*rot_shear_coef2*rot_shear_coef2*rot_shear_coef2);

    return {coord{trans_shear, 0, 0}, coord{0, 0, rot_shear}};

}

void acceleration(std::pair<coord, coord> Fbond, std::pair<coord, coord> Fshear, coord Frepulsion) {

    static double force_cal = 1 / (mass*0.001);
    static double force_cal2 = 1 / (rot_inertia*1e-12);
    np.acc = force_cal * (Fbond.first + Fshear.first + Frepulsion); // (N/kg)
    np.rot_acc = force_cal2 * (Fbond.second + Fshear.second + Frepulsion); // (s^-2)

}