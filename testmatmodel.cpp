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
    double x[12] = {3, 0, 0, 1, 0, 0, 3, 0, 0, 4, 0, 0};
    
    for(int i = -10; i < 10; i ++) {
        double h = mm.getLayerHeight(i, x);
        std::cout << i << ": " << h << "\n";
    }
    
    return 0;
}
