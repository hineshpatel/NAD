//
// Created by Mingqiu Wang on 5/12/16.
//
#include <iosfwd>
#include "declaration.h"

const coord operator+( const coord &a, const coord &b ) {
    return { a.x + b.x, a.y + b.y, a.z + b.z };
}
const coord operator-( const coord &a, const coord &b ) {
    return { a.x - b.x, a.y - b.y, a.z - b.z };
}
const coord operator*( const coord &a, double b ) {
    return coord{ a.x * b, a.y * b, a.z * b };
}
const coord operator*( double a, const coord &b ) {
    return coord{ b.x * a, b.y * a, b.z * a };
}
std::ostream& operator<<( std::ostream &os, const coord &vec) {
    os << vec.x << "\t" << vec.y << "\t" << vec.z;
    return os;
}