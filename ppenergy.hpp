/* 
 * File:   energy.hpp
 * Author: medved
 *
 * Energy computations
 * 
 * Created on August 24, 2015, 12:48 PM
 */

#ifndef PAIRPOTENTIALENERGY_HPP
#define	PAIRPOTENTIALENERGY_HPP

#include <math.h>
#include <functional>
#include "energy.hpp"
#include "matmodel.hpp"
#include <util/common/bnbdebug.hpp>

namespace lur {

    class PairPotentialEnergy : public Energy {
    public:

        /**
         * Constructor
         * @param model reference to the model
         * @param potent reference to the potential computing function         
         */
        PairPotentialEnergy(const MatModel& model, std::function <double (int, int, double) > potent) :
        mMatModel(model), mPotent(potent), mFixedAtoms(false), mOnceComputed(false) {
            mLBounds.resize(model.mNumLayers);
            mUBounds.resize(model.mNumLayers);
        }

        /**
         * Compute the energy of the material piece
         * @param x parameters, defining the layers height, displacements and interatomic distances
         * @return the energy
         */
        double energy(const double* x) {
            if (::bnbDebugVar0) {
                ::bnbDebugVar1 = 0;
            }
            double E = 0;
            int i;
#if 0            
            for (i = 0;; i++) {
                double h = mMatModel.getLayerHeight(i, x);
                if (h > mMatModel.mHeight)
                    break;
                E += layerEnergy(i, x);
            }
#endif            
            for (i = 0; i < mMatModel.mNumLayers; i++) {
                E += layerEnergy(i, x);
            }
            if (::bnbDebugVar0) {
                std::cout << "Contributed " << ::bnbDebugVar1 << " atoms from " << i << "layers \n";
            }
            return E;

        }

        /**
         * Retrieve math model
         * @return math model reference
         */
        const MatModel& getModel() const {
            return mMatModel;
        }

        /**
         * Sets fixed atoms regime on or off. Fixed atoms means that the energy is computed for a set of atoms rather than a piece of material
         * @param onoff on - true, off - false
         */
        void setFixedAtoms(bool onoff) {
            mFixedAtoms = onoff;
            mOnceComputed = false;
        }

    private:

        /**
         * Computes interaction energy for a layer
         * @param i layer's number
         * @param x layer's data
         * @return energy value
         */

        double layerEnergy(int i, const double* x) {
            double E = 0;
            if (!mFixedAtoms)
                mMatModel.computeBounds(x, mLBounds, mUBounds);
            else if (!mOnceComputed) {
                mOnceComputed = true;
                mMatModel.computeBounds(x, mLBounds, mUBounds);
            }
            for (int j = mLBounds[i]; j <= mUBounds[i]; j++) {
                E += atomEnergy(i, j, x);
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
        double atomEnergy(int i, int j, const double* x) const {
            double v = 0;
            ::bnbDebugVar1++;

            /**
             * Computes interaction energy for two atoms
             */
            auto inter = [&] (int i2, int j2) {

                if (!((i == i2) && (j == j2))) {
                    double q = mMatModel.getSqrDistance(i, j, i2, j2, x);
                    int a1 = mMatModel.mLayersAtoms[mMatModel.getReferenceLayer(i)];
                    int a2 = mMatModel.mLayersAtoms[mMatModel.getReferenceLayer(i2)];
                    v += mPotent(a1, a2, q);
                }
            };

            mMatModel.traverseLattice(i, j, inter, x);
            return v;
        }


        std::function< double (int, int, double) > mPotent;
        const MatModel& mMatModel;
        bool mFixedAtoms;
        bool mOnceComputed;
        std::vector< int > mLBounds;
        std::vector< int > mUBounds;

    };

}

#endif	/* ENERGY_HPP */

