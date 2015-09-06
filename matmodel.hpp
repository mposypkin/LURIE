/* 
 * File:   matmodel.hpp
 * Author: medved
 * 
 * Layered material model definition
 * 
 * Created on August 21, 2015, 1:30 PM
 */

#ifndef MATMODEL_HPP
#define	MATMODEL_HPP

#include <math.h>

namespace lur {

    static const unsigned int maxLayers = 100;

    /**
     * Material model definition
     */
    struct MatModel {

        /**
         * Computes layers height by its number
         * @param i the number (may be greater or less than zero) 
         * @param x the layer's data
         * @return layer's height from the 0 level (may be positive, zero or negative)
         */
        double getLayerHeight(int i, const double* x) {
            double rv;
            if (i == 0)
                rv = 0;
            else {
                auto getla = [x](int I) {
                    return (I == 0) ? 0 : x[3 * I];
                };
                int k = mNumLayers;
                double h = getla(k - 1) + x[0];
                int m = floor((double) i / (double) k);
                int r = i - m * k;
                rv = getla(r) + h * m;
            }
            
            return rv;
        }
        
        /**
         * Number of main layers
         */
        int mNumLayers;

        /**
         * The maximal radius of potential interaction
         */
        double mRadius;

        /**
         * The length of the material piece to model
         */
        double mLength;

        /**
         * Layer's atoms, each atom is identified by an integral number
         */
        int mLayersAtoms[maxLayers];


    };
}

#endif	/* MATMODEL_HPP */

