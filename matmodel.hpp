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
#include <iostream>

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
        double getLayerHeight(int i, const double* x) const {
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
         * Computes the reference layer in model data for an arbitrary row
         * taking into account periodicity
         * @param i row number
         * @return reference layer
         */
        int getReferenceLayer(int i) const {
            int k = mNumLayers;
            int m = floor((double) i / (double) k);
            return i - m * k;
        }
        
        /**
         * Computes the displacement of the first atom and the interatomic distance in a layer
         * @param l layers number
         * @param d displacment
         * @param x layers data
         * @param s stride
         */
        void getDisplacementAndStride(int l, double& d, double& s, const double* x) const {
            int r = getReferenceLayer(l);            
            d = x[3 * r + 1];
            s = x[3 * r + 2];
        }
        
        /**
         * Retrieve atoms 'x' coordinate (positive or negative offset from zero)
         * @param i atoms layer
         * @param j atoms number
         * @param x layers data
         * @return offset
         */
        double getOffset(int i, int j, const double* x) const {
            double d, s;
            getDisplacementAndStride(i, d, s, x);
            return d + j * s;
        }

        
        /**
         * Computes the square of euclidian distance between two atoms given by their 
         * layer number and number in a layer
         * @param i1 first atoms layer
         * @param j1 first atoms number
         * @param i2 second atoms layer
         * @param j2 second atoms number
         * @return distance
         */
        double getSqrDistance(int i1, int j1, int i2, int j2, const double* x) const {
            double y1 = getLayerHeight(i1, x);
            double x1 = getOffset(i1, j1, x);
            double y2 = getLayerHeight(i2, x);
            double x2 = getOffset(i2, j2, x);
            return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
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
         * The height of the material piece to model
         */
        double mHeight;
        
        /**
         * Layer's atoms, each atom is identified by an integral number
         */
        int mLayersAtoms[maxLayers];


    };
}

#endif	/* MATMODEL_HPP */

