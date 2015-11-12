/* 
 * File:   energy.hpp
 * Author: medved
 *
 * Created on October 23, 2015, 5:01 PM
 */

#ifndef ENERGY_HPP
#define	ENERGY_HPP

namespace lur {

    /**
     * Abstract class of the energy of the material piece
     */
    class Energy {
    public:

        /**
         * Compute the energy of the material piece
         * @param x parameters, defining the layers height, displacements and interatomic distances
         * @return the energy
         */
        virtual double energy(const double* x) = 0;
    };
}

#endif	/* ENERGY_HPP */

