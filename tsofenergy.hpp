/* 
 * File:   tsofenergy.hpp
 * Author: medved
 *
 * Created on October 24, 2015, 11:59 AM
 */

#ifndef TSOFENERGY_HPP
#define	TSOFENERGY_HPP

#include <math.h>
#include <functional>
#include <vector>
#include <util/common/bnbdebug.hpp>
#include "energy.hpp"
#include "matmodel.hpp"
#include "tersoffparams.hpp"
#include "tersoffutils.hpp"

namespace lur {

    /**
     * Lattice energy by tersoff potential
     */
    class TersoffEnergy : public Energy {
    public:

        /**
         * Constructor
         * @param model reference to the model   
         */
        TersoffEnergy(const MatModel& model, const TersoffUtils& tutils) :
        mMatModel(model), mTutils(tutils), mFixedAtoms(false), mOnceComputed(false) {
            mLBounds.resize(model.mNumLayers);
            mUBounds.resize(model.mNumLayers);
        }

        /**
         * Compute the energy of the material piece
         * @param x parameters, defining the layers height, displacements and interatomic distances
         * @return the energy
         */
        double energy(const double* x)  {
            double E = 0;
#if 0            
            for (int i = 0;; i++) {
                double h = mMatModel.getLayerHeight(i, x);
                if (h > mMatModel.mHeight)
                    break;
                E += layerEnergy(i, x);
            }
#endif            
            for (int i = 0; i < mMatModel.mNumLayers; i++) {
                E += layerEnergy(i, x);
            }

            return E;
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
        double layerEnergy(int i, const double* x)  {
            double E = 0;
            if (!mFixedAtoms)
                mMatModel.computeBounds(x, mLBounds, mUBounds);
            else if(!mOnceComputed) {
                mOnceComputed = true;
                mMatModel.computeBounds(x, mLBounds, mUBounds);
            }
            for (int j = mLBounds[i]; j <= mUBounds[i]; j++) {
                E += atomEnergy(i, j, x);
            }
            return E;
        }

        /**
         * Computes interaction energy for an atom 
         * @param i atoms' layer
         * @param j atoms' number in a layer
         * @param x layer's data
         * @return energy value
         */
        double atomEnergy(int i, int j, const double* x) const {
            double v = 0;


            /**
             * Computes interaction energy for two atoms
             */
            auto inter = [&] (int i2, int j2) {
                if (!((i == i2) && (j == j2))) {
                    double q = mMatModel.getSqrDistance(i, j, i2, j2, x);
                    double r = sqrt(q);
                    double fc = mTutils.cutoff(r);
                    if (fc > 0) {
                        double zij = 0;

                        /**
                         * Aux function to compute angular term zij
                         */
                        auto ang = [&] (int i3, int j3) {
                            //std::cout << "i3 = " << i3 << "j3= " << j3 << "\n";
                            bool isij = (i == i3) && (j == j3);
                            bool isi2j2 = (i2 == i3) && (j2 == j3);
                            if (!(isij || isi2j2)) {
                                double q3 = mMatModel.getSqrDistance(i, j, i3, j3, x);
                                double r3 = sqrt(q3);
                                double fc3 = mTutils.cutoff(r3);
                                if (fc3 > 0) {
                                    double o = mTutils.computeOmega(r, r3);
                                    double scalv = mMatModel.getScalarMult(i, j, i2, j2, i3, j3, x);
                                    double cosv = scalv / (r * r3);
                                    double u = mTutils.computeG(cosv);
                                    zij += fc3 * u * o;
                                }
                            }
                        };
                        mMatModel.traverseLattice(i, j, ang, x);
                        double bij = mTutils.computeBij(zij);
                        v += fc * (mTutils.VR(r) - bij * mTutils.VA(r));
                    }
                }
            };

            mMatModel.traverseLattice(i, j, inter, x);
            return v;
        }

        const TersoffUtils& mTutils;
        const MatModel& mMatModel;
        bool mFixedAtoms;
        bool mOnceComputed;

        std::vector< int > mLBounds;
        std::vector< int > mUBounds;
    };

}


#endif	/* TSOFENERGY_HPP */

