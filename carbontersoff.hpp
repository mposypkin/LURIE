/* 
 * File:   carbontersoff.hpp
 * Author: medved
 *
 * Created on October 24, 2015, 4:59 PM
 */

#ifndef CARBONTERSOFF_HPP
#define	CARBONTERSOFF_HPP

#include "tersoffparams.hpp"

namespace lur {

    /**
     * Based on Powell D. Elasticity, Lattice Dynamics and Parameterisation 
     * Techniques for the Tersoff Potential Applied to Elemental and 
     * Type III-V Semiconductors :University of Sheffield, 2006.
     * @param tparam Tersoff structure
     */
    static void fillCarbonParametersPowell(TersoffParams& tparam) {
        tparam.mDe = 5.164533;
        tparam.mS = 1.576879;
        tparam.mB = 1.964037;
        tparam.mRe = 1.447206;
        tparam.mR = 1.95;
        tparam.mRcut = 0.15;
        tparam.mC = 38049;
        tparam.mD = 4.384;
        tparam.mH = -0.5758;
        tparam.mNu = 0.72751;
        tparam.mGamma = 1.5724E-7;
        tparam.mLambda = 1.5;
    }

    /**
     * Original Tersoff parameters from 
     * Tersoff J. Empirical interatomic potential for carbon, with applications 
     * to amorphous carbon Physical Review Letters. 1988. V. 61. 25. 
     * @param tparam Tersoff structure
     */
    static void fillCarbonParametersTersoffOriginal(TersoffParams& tparam) {
        tparam.mDe = 5.166159;
        tparam.mS = 1.576879;
        tparam.mB = 1.964037;
        tparam.mRe = 1.4471159;
        tparam.mR = 1.95;
        tparam.mRcut = 0.15;
        tparam.mC = 38049;
        tparam.mD = 4.3484;
        tparam.mH = -0.57058;
        tparam.mNu = 0.72751;
        tparam.mGamma = 1.5724E-7;
        tparam.mLambda = 0;
    }

    /**
     * Modified Tersoff parameters from 
     * Lindsay, L., & Broido, D. A. (2010). Optimized Tersoff and Brenner empirical potential parameters for lattice 
     * dynamics and phonon thermal transport in carbon nanotubes and graphene. Physical Review B, 81(20), 205441. 
     * @param tparam Tersoff structure
     */
    static void fillCarbonParametersTersoffOptimized(TersoffParams& tparam) {
        tparam.mDe = 9.303600;
        tparam.mS = 1.576879;
        tparam.mB = 1.964037;
        tparam.mRe = 1.278456;
        tparam.mR = 1.95;
        tparam.mRcut = 0.15;
        tparam.mC = 38049;
        tparam.mD = 4.3484;
        tparam.mH = -0.930;
        tparam.mNu = 0.72751;
        tparam.mGamma = 1.5724E-7;
        tparam.mLambda = 0;
    }

}

#endif	/* CARBONTERSOFF_HPP */

