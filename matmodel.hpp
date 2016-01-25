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
#include <functional>
#include <util/common/vec.hpp>
#include <util/common/bnbdebug.hpp>

namespace lur {

    static const unsigned int maxLayers = 100;

    /**
     * Material model definition
     */
    struct MatModel {

        /**
         * Computes layers height by its number
         * @param i the number of the layer
         * @param x the layer's data
         * @return layer's height from the 0 level (may be positive, zero or negative)
         */
        double getLayerHeight(int i, const double* x) const {
            double h = 0;
            if (i > 0) {
                for (int j = 1; j <= i; j++) {
                    int J = getReferenceLayer(j);
                    h += x[3 * J];
                }
            } else if (i < 0) {
                for (int j = 0; j > i; j--) {
                    int J = getReferenceLayer(j);
                    h -= x[3 * J];
                }
            }
            return h;
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
         * Computes scalar multiple of vectors r12 and r13
         * @param i1 first atoms layer
         * @param j1 first atoms number
         * @param i2 second atoms layer
         * @param j2 second atoms number
         * @param i3 third atoms layer
         * @param j3 third atoms number
         * @param x
         * @return 
         */
        double getScalarMult(int i1, int j1, int i2, int j2, int i3, int j3, const double* x) const {
            double y1 = getLayerHeight(i1, x);
            double x1 = getOffset(i1, j1, x);
            double y2 = getLayerHeight(i2, x);
            double x2 = getOffset(i2, j2, x);
            double y3 = getLayerHeight(i3, x);
            double x3 = getOffset(i3, j3, x);
            double scalv = (x2 - x1) * (x3 - x1) + (y2 - y1) * (y3 - y1);
            return scalv;
        }

        /**
         * Traverses the neighbourhood (including atom itself) of a given atom defined by the model cut radius
         * @param i atom's layer
         * @param j atom's position in a layer
         * @param compf functor to call for each atom
         */
        void traverseLattice(int i, int j, std::function <void (int, int) > compf, const double* x) const {
            double myx = getOffset(i, j, x);
            double myy = getLayerHeight(i, x);
            /**
             * Computes interaction energy with the layer
             */
            auto lenerg = [&] (int l) {
                double d;
                double s;
                getDisplacementAndStride(l, d, s, x);
                int tl = ceil((myx - mRadius - d) / s);
                int tu = floor((myx + mRadius - d) / s);
                double u = 0;
                for (int t = tl; t <= tu; t++) {
                    compf(l, t);
                }
                return u;
            };


            /**
             * Add energy of 'my' layer
             */
            lenerg(i);
            /**
             * Pass up              
             */
            for (int l = i + 1;; l++) {
                double h = getLayerHeight(l, x);
                if (h - myy > mRadius)
                    break;
                lenerg(l);
            }
            /**
             * Add pass down
             */
            for (int l = i - 1;; l--) {
                double h = getLayerHeight(l, x);
                if (myy - h > mRadius)
                    break;
                lenerg(l);
            }
        };

        /**
         * Computes bounds of atom indeces based on the rectangle shape
         * @param x lattice parameters
         * @param lbounds lower indices
         * @param rbounds upper indices
         */
        void computeBounds(const double* x, std::vector<int> &lbounds, std::vector<int> &rbounds) const {
            for (int i = 0; i < mNumLayers; i++) {
                int a = 0, b = -1;
                int j = -1;
                while (getOffset(i, j, x) >= 0) {
                    a = j;
                    j--;
                }
                j = 0;
                while (getOffset(i, j, x) <= mLength) {
                    b = j;
                    j++;
                }
                lbounds[i] = a;
                rbounds[i] = b;
            }
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

