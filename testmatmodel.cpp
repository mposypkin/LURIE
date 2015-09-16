/* 
 * File:   testmatmodel.cpp
 * Author: medved
 *
 * Created on September 6, 2015, 2:36 PM
 */

#include <iostream>
#include "matmodel.hpp"

/*
 * Testing mathematical model
 */
int main(int argc, char** argv) {

    lur::MatModel mm;
    
    mm.mNumLayers = 4;
    double x[12] = {3, 1, 1,   1, 2, 1,   3, 1, 1,   4, 2, 1};
    
    for(int i = -10; i < 10; i ++) {
        double h = mm.getLayerHeight(i, x);
        std::cout << i << ": " << h << "\n";
    }
    
    double d = mm.getSqrDistance(1, 4, 5, 4, x);
    std::cout << " d = " << d << "\n";
    return 0;
}
