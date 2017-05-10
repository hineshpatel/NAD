//
// Created by Mingqiu Wang on 5/7/16.
//

#ifndef NANOAD_DATASTRUCTURES_H
#define NANOAD_DATASTRUCTURES_H
#include <ostream>

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

    coord position;
    bool bound;
    int pair;

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

    coord position;
    coord position_origin;
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

    cutoff() {}
    cutoff(double bondL, double ratio) : deltaMax{ratio * bondL}, bondLMax{ratio * bondL + bondL},
                                         bondLMin{ratio * bondL + bondL} {}
};

#endif //NANOAD_DATASTRUCTURES_H
