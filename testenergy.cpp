/* 
 * File:   testenergy.cpp
 * Author: medved
 *
 * Created on August 21, 2015, 2:03 PM
 */
#include <assert.h>
#include <functional>
#include <util/common/bnberrcheck.hpp>

#include "matmodel.hpp"
#include "ppenergy.hpp"
#include "parsejson.hpp"

/*
 * 
 */

static char models [] = "{"
    "\"model\": {"
        "\"nlay\": 4,"
        "\"range\": 1.9,"
        "\"length\": 5.5,"
        "\"height\": 6.5,"
        "\"atoms\": [1, 1, 1, 1]"
   "},"
   "\"lattice\": {\"v\": 0, \"x\": [1, 0, 2,  1, 1, 2,  1, 0, 2,  1, 1, 2]}"
"}";

double sqrpotent(int a1, int a2, double q) {
    return (q - 1)*(q - 1);
}

int main(int argc, char** argv) {
    lur::MatModel mm;
    lur::ParseJson::parseModelData(models, mm);
    /*
    mm.mNumLayers = 4;
    mm.mLength = 5.5;
    mm.mHeight = 6.5;
    mm.mRadius = 1.9;
     */
    std::cout << "mm.mNumLayers = " << mm.mNumLayers << "\n";
    std::cout << "mm.mLength = " << mm.mLength << "\n";
    std::cout << "mm.mHeight = " << mm.mHeight << "\n";
    std::cout << "mm.mRadius = " << mm.mRadius << "\n";
    lur::PairPotentialEnergy enrg(mm, sqrpotent);
    //double x[12] = {1, 0, 2, 1, 1, 2, 2, 0, 2, 3, 1, 2};
    double x[12];
    double v;
    lur::ParseJson::parseLatticeData(models, mm, v, x);
    double E = enrg.energy(x);
    std::cout << E << "\n";
    BNB_ASSERT(E == 84);
    mm.mRadius = 2.9;
    double y[12] = {2, 0, 2,  2, 1, 2,  2, 0, 2,  2, 1, 2};
    E = enrg.energy(y);
    std::cout << E << "\n";
    BNB_ASSERT(E == 984);
    return 0;
}

