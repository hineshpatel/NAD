////
//// Created by Mingqiu Wang on 6/27/16.
////
//
//#include "declaration.h"
//#include "dataStructures.h"
//using namespace std;
//
//
///**
// * Test forces and acceleration
// *
// * @func: Fbond(), Fshear(), acceleration(), Frepulsion()
// *
// * @param
// * @return: acc and rot_acc
// *
// */
//pair<coord, coord> force_and_acc( int bond_num, double par[],
//                                  double rec[][3], double lig[][3]) {
//
//
//    np.position = {par[0], par[1], par[2]};
//    for (auto x = 0; x < bond_num; ++x) {
//        activeBonds.insert(x);
//        receptor receptor1;
//        receptor1.position = {rec[x][0], rec[x][1], rec[x][2]};
//        receptors.push_back(receptor1);
//        ligand ligand1;
//        ligand1.updatePA({lig[x][0],lig[x][1],lig[x][2]}, np.position);
//        ligands.push_back(ligand1);
//        bond bond1;
//        bond1.ligand = x;
//        bond1.receptor = x;
//        bond1.delta = dist(ligands.at(x).position, receptors.at(x).position) - _bondEL;
//        bonds.push_back(bond1);
//    }
//
//
//    std::pair<coord, coord> bond = Fbond(activeBonds, bonds, receptors, ligands);
//    std::pair<coord, coord> shear = Fshear(np.position.z);
//    acceleration(bond, shear, Frepulsion(np.position.z));
//    return {np.acc, np.rot_acc};
//}
//
///**
// * Test translation and rotation
// *
// * @func: translation(), rotateLig(), rotate()
// *
// * @param
// * @return:
// *
// */
//void trans_and_rot( double v[], double rot_v[]) {
//
//    setRNG(sfmt);
//
//    np.velocity = {v[0], v[1], v[2]};
//    np.rot_velocity = {rot_v[0], rot_v[1], rot_v[2]};
//    translation(np.velocity, np.position, np.acc); // translate nanoparticle
//    rotateLig(ligands, rotate(np.rot_velocity, np.rot_acc), np.position); // rotate nanoparticle
//
//}
//
//
//
//
//
//int main() {
//
//// test case 1 ================================================================================
////    int bond_num = 3;
////    double par[] = {1.553678, -45.818843, 124.235014};
////    double rec[][3] = {0, -45.0, 0, 3, -50, 0, -0.5, -40, 0};
////    double lig[][3] = {58.055188,-49.81867,35.823587, 77.100255,15.987777,85.535451,
////                       -9.222984,-149.624129,135.781703};
//// end input ================================================================================
//
////// test case 2 ==============================================================================
////
////    int bond_num = 5;
////    double par[] = {226.475992,-132.817602,171.916843};
////    double rec[][3] = {-1014.335742, -672.660289, 0,
////                       -818.888729, 233.272184, 0,
////                       759.950703, 757.620378, 0,
////                       560.392341, -54.561730, 0,
////                       -682.475808, 216.211236, 0};
////    double lig[][3] = {316.135123,	-86.091697,	143.580827,
////                       247.240758,	-46.888818,	228.574309,
////                       250.868839,	-49.537881,	231.030944,
////                       207.489674,	-97.504291,	268.960595,
////                       210.489308,	-104.408229,	271.728332};
////
////// end input ================================================================================
//
////    force_and_acc(bond_num, par, rec, lig);
////    cout << "particle trans and rot accelerations\t" << np.acc << "\t" << np.rot_acc << endl;
////
//
//// test case 1 =============================================================================
////    double v[] = {28533362.197028,23205777.429526,-39590677.153833};
////    double rot_v[] = {478586.526767,-136219.076178,603880.867970};
//// end input ================================================================================
//
////    trans_and_rot(v, rot_v);
////    printf("particle and lig pos %e\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", np.position.x, np.position.y,
////               np.position.z, np.velocity.x, ligands.at(0).position.x, ligands.at(1).position.x, ligands.at(2).position.x);
////
////
////
////// test seg-seg distance
////    cout << seg_seg_Dist(
////            {0.000000e+00, 0.000000e+00, 0.000000e+00},
////            {2.278302e+01, 1.695014e+01, 2.975076e+01},
////            {-1.580203e+00, 3.034802e-01, 0.000000e+00},
////            {2.359700e+01, 1.263818e+01, 3.020311e+01}
////    ) << endl;
//
//
//// test putNPBack
//    coord NPposition = {-16, 0, 5};
//    putNPBack(NPposition);
//    cout << NPposition << endl;
//}
//
//
//
