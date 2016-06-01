//
// Created by Mingqiu Wang on 5/7/16.
//

#ifndef NANOAD_DATASTRUCTURES_H
#define NANOAD_DATASTRUCTURES_H


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


class receptor {

public:

    coord position;
    bool bind;
    int pair;

    receptor() : bind {false}, pair {-1} {}

    void unpairing() {
        bind = false;
        pair = -1;
    }

    void pairing(int i) {
        bind = true;
        pair = i;
    }
};

class ligand {

public:

    coord position;
    coord position_origin;
    bool bind;
    int pair;

    ligand() : bind {false}, pair {-1} {}

    void updatePO(const coord& origin, const coord& np) {
        position_origin = origin;
        position = position_origin + np;
    }

    void updatePA(const coord& pos, const coord& np) {
        position = pos;
        position_origin = position - np;
    }
    void unpairing() {
        bind = false;
        pair = -1;
    }

    void pairing(int i) {
        bind = true;
        pair = i;
    }

};

struct np {
    int name = 0;
    coord position;
    coord lastPairPos;
    coord angle;
    coord velocity;
    coord rot_velocity;
    coord acc;
    coord rot_acc;
};

struct bond {

    int name = 0;
    bool bind = false;
    int ligand = 0;
    int receptor = 0;

    double formTime = 0; // (s)
    coord formPositionLigand; // (nm)
    coord formPositionReceptor; // (nm)

    double breakTime = 0; // (s)
    coord breakPositionLigand; // (nm)
    coord breakPositionReceptor; // (nm)

};

struct cutoff {

    double bondLMax;
    double bondLMin;
    double deltaMax;

};

#endif //NANOAD_DATASTRUCTURES_H
