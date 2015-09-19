/* 
 * File:   testenergy.cpp
 * Author: medved
 *
 * Created on August 21, 2015, 2:03 PM
 */
#include <assert.h>
#include <functional>

#include "matmodel.hpp"
#include "energy.hpp"

/*
 * 
 */


double sqrpotent(int a1, int a2, double q) {
    return (q - 1)*(q - 1);
}

int main(int argc, char** argv) {
    lur::MatModel mm;
    mm.mNumLayers = 4;
    mm.mLength = 5.5;
    mm.mHeight = 6.5;
    mm.mRadius = 1.9;
    lur::Energy enrg(mm, sqrpotent);
    double x[12] = {1, 0, 2, 1, 1, 2, 2, 0, 2, 3, 1, 2};
    double E = enrg.energy(x);
    std::cout << E << "\n";
    assert(E == 84);
    mm.mRadius = 2.9;
    double y[12] = {2, 0, 2, 2, 1, 2, 4, 0, 2, 6, 1, 2};
    E = enrg.energy(y);
    std::cout << E << "\n";
    assert(E == 984);
    return 0;
}

