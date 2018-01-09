//
// Created by Mingqiu Wang on 5/7/16.
//

#ifndef NANOAD_DATASTRUCTURES_H
#define NANOAD_DATASTRUCTURES_H
#include <ostream>
#include "parameters.h"

struct coord {
    double x;
    double y;
    double z;
    coord() : x{0}, y{0}, z{0} {}
    coord(double xx, double yy, double zz): x{xx}, y{yy}, z{zz} {}
};

const coord operator+( const coord &a, const coord &b );
const coord operator-( const coord &a, const coord &b );
const coord operator*( const coord &a, double b );
const coord operator*( double a, const coord &b );
std::ostream& operator<<( std::ostream &os, const coord &vec);


class receptor {

public:

    coord position; // if ORI is enabled, this is the position of tip of ICAM-1
    // otherwise it's the position of stem of ICAM-1 on surface
    bool bound;
    int pair;
    coord stem; // stem of ICAM-1 on surface

    receptor() : bound{false}, pair {-1} {}

    void unpairing() {
        bound = false;
        pair = -1;
    }

    void pairing(int i) {
        bound = true;
        pair = i;
    }
};

class ligand {

public:

    coord position; // if ORI is enabled, this is the position of actual Fab tip,
    // otherwise it's the actual position of Ab stem on the NP
    coord position_origin; // if ORI is enabled, this is the position of origin-based Fab tip,
    // otherwise it's the origin-based position of Ab stem on the NP
    bool bound;
    int pair;

    ligand() : bound{false}, pair {-1} {}

    void updatePO(const coord& origin, const coord& np) {
        position_origin = origin;
        position = position_origin + np;
    }

    void updatePA(const coord& pos, const coord& np) {
        position = pos;
        position_origin = position - np;
    }
    void unpairing() {
        bound = false;
        pair = -1;
    }

    void pairing(int i) {
        bound = true;
        pair = i;
    }

};

struct np {
    int name = 0;
    coord position;
    coord lastPairPos;
    coord velocity;
    coord rot_velocity;
    coord acc;
    coord rot_acc;
    coord iniPos;
    double start_time = 0.0;
};

struct bond {

    int name = -1;
    bool bound = false;
    int ligand = -1;
    int receptor = -1;

    double formTime = 0; // (s)
    coord formPositionLigand; // (nm)
    coord formPositionReceptor; // (nm)

    double breakTime = -1; // (s)
    coord breakPositionLigand; // (nm)
    coord breakPositionReceptor; // (nm)

    double delta = -1.0; // (nm) = current bond length - equilibrium length; (+) extension (-) compression

};

struct cutoff {

    double bondLMax; // max bond length (nm)
    double bondLMin; // min bond length (nm)
    double deltaMax; // max bond length - equilibrium bond length (nm)
    double fabMax; // max fab length (nm)
    double fabMin; // min fab length (nm)
    cutoff() {}
    cutoff(double bondL, double ratio, double fab_length = 6.4) :
            deltaMax{ratio * bondL},
            bondLMax{ratio * bondL + bondL},
            bondLMin{-1.0 * ratio * bondL + bondL},
            fabMax{ratio * bondL + fab_length},
            fabMin{-1.0 * ratio * bondL + fab_length}
    {}
};

#endif //NANOAD_DATASTRUCTURES_H
