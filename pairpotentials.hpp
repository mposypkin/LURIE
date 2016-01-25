/**
 * File:   pairpotentials.hpp
 * Author: medved
 *
 * Created on January 6, 2016, 5:48 PM
 * 
 * Pair potential functions that could be used to evaluate the energy of material piece
 */

#ifndef PAIRPOTENTIALS_HPP
#define	PAIRPOTENTIALS_HPP

namespace lur {

    /**
     * Square potential
     * @param a1 the type of the first atom
     * @param a2 the type of the second atom
     * @param q the square distance
     * @return the value of potential
     */
    double sqrpotent(int a1, int a2, double q) {
        return (q - 1)*(q - 1) - 1;
    }

    /**
     * Lennard-Jones potential
     * @param a1 the type of the first atom
     * @param a2 the type of the second atom
     * @param q the square distance
     * @return the value of potential
     */
    double ljpotent(int a1, int a2, double q) {
        if (q == 0)
            q = 0.00001;
        double u = q * q * q;
        double v = u * u;
        double p = 1. / v - 2. / u;
        return p;
    }

    /**
     * Morse potential
     * @param a1 the type of the first atom
     * @param a2 the type of the second atom
     * @param q the square distance
     * @return the value of potential
     */
    double morsepotent(int a1, int a2, double q) {
        BNB_ASSERT(q >= 0);
        double r = sqrt(q);
        const double rho = 14;
        double E = exp(rho * (1 - r));
        double p = E * (E - 2);
        return p;
    }

    /**
     * Cuts off the pair potential. Delivered smoothely (up to 2nd derivative)
     * cut potential from Q-D to Q
     */
    struct PotentialCutter {

        /**
         * Constructor
         * @param Q cutoff square radius
         * @param D cutoff delta
         * @param potent potential function
         */
        PotentialCutter(double Q, double D, std::function<double (int, int, double)> potent)
        : mQ(Q), mD(D), mPotent(potent) {
        }

        /**
         * Truncated potential
         * @param a1 the type of the first atom
         * @param a2 the type of the second atom
         * @param q the square distance
         * @return the value of potential
         */
        double operator () (int a1, int a2, double q) {
            double m;
            if (q <= mQ - mD) {
                return mPotent(a1, a2, q);
            } else if (q >= mQ) {
                return 0;
            } else {
                auto z = [](double x) {
                    double x3 = x * x * x;
                    double x5 = x3 * x * x;
                    return -3. / 16. * x5 + 5. / 8. * x3 - 15. / 16. * x + 0.5;
                };
                double y = (2. / mD) * q + 1 - (2. * mQ / mD);
                m = z(y);
                return m * mPotent(a1, a2, q);
            }
        }

        double mQ;
        double mD;
        std::function<double (int a1, int a2, double q) > mPotent;

    };
}

#endif	/* PAIRPOTENTIALS_HPP */

