/* 
 * File:   tersoffutils.hpp
 * Author: medved
 *
 * Created on October 24, 2015, 10:04 AM
 */

#ifndef TERSOFFUTILS_HPP
#define	TERSOFFUTILS_HPP

#include <math.h>
#include "tersoffparams.hpp"

namespace lur {

    /**
     * Means for working with Tersoff potential
     */
    class TersoffUtils {
    public:

        /**
         * Constructor
         * @param tersp reference to tersoff parameters
         */
        TersoffUtils(const TersoffParams& tersp) : mTpar(tersp) {
        }

        /**
         * Cutoff function
         * @param r radius
         * @return cut value between 0 and 1
         */
        double cutoff(double r) const {
            double rv;
            if (r < mTpar.mR - mTpar.mRcut) {
                rv = 1;
            } else if (r > mTpar.mR + mTpar.mRcut) {
                rv = 0;
            } else {
                rv = 0.5 * (1 - sin(M_PI * (r - mTpar.mR) / (2. * mTpar.mRcut)));
            }
            return rv;
        }

        /**
         * Repulsive energy
         * @param r distance
         * @return repulsive energy
         */
        double VR(double r) const {
            double rv;
            rv = mTpar.mDe / (mTpar.mS - 1.);
            rv *= exp(-mTpar.mB * sqrt(2. * mTpar.mS) * (r - mTpar.mRe));
            return rv;
        }

        /**
         * Attractive energy
         * @param r distance
         * @return attractive energy
         */
        double VA(double r) const {
            double rv;
            rv = mTpar.mDe * mTpar.mS / (mTpar.mS - 1.);
            rv *= exp(-mTpar.mB * sqrt(2. / mTpar.mS) * (r - mTpar.mRe));
            return rv;
        }

        /**
         * Computes angular term
         * @param costheta the cosinus of the angle
         */
        double computeG(double costheta) const {
            double rv = 1;
            double v = mTpar.mC / mTpar.mD;
            rv += v * v;
            v = mTpar.mH - costheta;
            rv -= mTpar.mC * mTpar.mC / (mTpar.mD * mTpar.mD + v * v);
            return rv;
        }

        /**
         * Compute bij for i,j term
         * @param zeta zeta for i, j term
         * @return value of bij
         */
        double computeBij(double zeta) const {
            double rv;
            double a = 1. + pow(mTpar.mGamma * zeta, mTpar.mNu);
            rv = pow(a, -1. / (2. * mTpar.mNu));
            return rv;
        }

        /**
         * Computes Omega ik term
         * @param rij distance between  atoms i and j
         * @param rik distance between atoms i and k
         * @return value of omega
         */
        double computeOmega(double rij, double rik) const {
            double v = mTpar.mLambda * (rij - rik);
            double w = exp(v * v * v);
            return w;
        }
        
    private:
        const TersoffParams& mTpar;
    };
}

#endif	/* TERSOFFUTILS_HPP */

