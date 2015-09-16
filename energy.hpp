/* 
 * File:   energy.hpp
 * Author: medved
 *
 * Energy computations
 * 
 * Created on August 24, 2015, 12:48 PM
 */

#ifndef ENERGY_HPP
#define	ENERGY_HPP

#include <math.h>
#include <functional>
#include "matmodel.hpp"

namespace lur {

    class Energy {
    public:

        /**
         * Constructor
         * @param model reference to the model
         * @param potent reference to the potential computing function         
         */
        Energy(const MatModel& model, std::function <double (int, int, double) > potent) :
        mMatModel(model), mPotent(potent) {
        }

        /**
         * Compute the energy of the material piece
         * @param x parameters, defining the layers height, displacements and interatomic distances
         * @return the energy
         */
        double energy(const double* x) {
            double E = 0;
            /* Atoms are numbered from 0 */
            for (int i = 0;; i++) {
                double h = mMatModel.getLayerHeight(i, x);
                if (h > mMatModel.mHeight)
                    break;
                E += layerEnergy(i, x);
            }
            return E;

        }

    private:

        /**
         * Computes interaction energy for a layer
         * @param i layer's number
         * @param x layer's data
         * @return energy value
         */
        double layerEnergy(int i, const double* x) {
            int j = 0;
            double E = 0;
            while (mMatModel.getOffset(i, j, x) < mMatModel.mLength) {
                E += atomEnergy(i, j, x);
                j++;
            }
            return E;
        }

        /**
         * Computes interaction energy for an atom (interactions with all neighbours are accounted)
         * @param i atoms' layer
         * @param j atoms' number in a layer
         * @param x layer's data
         * @return energy value
         */
        double atomEnergy(int i, int j, const double* x) {
            double y = mMatModel.getLayerHeight(i, x);
            double v = 0;
            double R = mMatModel.mRadius;
            double myx = mMatModel.getOffset(i, j, x);

            /**
             * Computes interaction energy for two atoms
             */
            auto inter = [&] (int i1, int j1, int i2, int j2) {
                double q = mMatModel.getSqrDistance(i1, j1, i2, j2, x);
                int a1 = mMatModel.mLayersAtoms[mMatModel.getReferenceLayer(i1)];
                int a2 = mMatModel.mLayersAtoms[mMatModel.getReferenceLayer(i2)];
                double v = mPotent(a1, a2, q);
                return v;
            };

            /**
             * Computes interaction energy with the layer
             */
            auto lenerg = [&] (int l) {
                double d;
                double s;
                mMatModel.getDisplacementAndStride(l, d, s, x);
                int tl = ceil((myx - R - d) / s);
                int tu = floor((myx + R - d) / s);
                double u = 0;
                if (i != l) {
                    for (int t = tl; t <= tu; t++) {
                        u += inter(i, j, l, t);
                    }
                } else {
                    for (int t = tl; t <= tu; t++) {
                        if (j != t)
                            u += inter(i, j, l, t);
                    }
                }
                return u;
            };
            /**
             * Add energy of 'my' layer
             */
            v = lenerg(i);
            /**
             * Pass up              
             */
            for (int l = i + 1;; l++) {
                double h = mMatModel.getLayerHeight(l, x);
                if (h - y > R)
                    break;
                v += lenerg(l);
            }
            /**
             * Add pass down
             */
            for (int l = i - 1;; l--) {
                double h = mMatModel.getLayerHeight(l, x);
                if (y - h > R)
                    break;
                v += lenerg(l);
            }

            return v;
        }


        std::function< double (int, int, double) > mPotent;
        const MatModel& mMatModel;
    };

}

#endif	/* ENERGY_HPP */

