/* 
 * File:   tersoffparams.hpp
 * Author: medved
 *
 * Created on October 24, 2015, 10:01 AM
 */

#ifndef TERSOFFPARAMS_HPP
#define	TERSOFFPARAMS_HPP

namespace lur {

    /**
     * Tersoff parameters
     */
    struct TersoffParams {
        /**
         * Pair interaction energy
         */
        double mDe;

        /**
         * Pair interaction separation
         */
        double mRe;

        /**
         * Beta (pair interaction)
         */
        double mB;

        /**
         * S (pair interaction)
         */
        double mS;

        /**
         * n parameter
         */
        double mNu;

        /**
         * gamma parameter
         */
        double mGamma;

        /**
         * Lambda parameter
         */
        double mLambda;

        /**
         * c parameter
         */
        double mC;

        /**
         * d parameter
         */
        double mD;

        /**
         * h parameter
         */
        double mH;

        /**
         * Cutoff center
         */
        double mR;

        /**
         * Cutoff radius
         */
        double mRcut;

    };

}

#endif	/* TERSOFFPARAMS_HPP */

