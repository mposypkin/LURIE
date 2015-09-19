/* 
 * File:   lurobj.hpp
 * Author: medved
 *
 * Objective function for material model
 * 
 * Created on September 18, 2015, 10:32 AM
 */

#ifndef LUROBJ_HPP
#define	LUROBJ_HPP

#include <problems/optlib/objective.hpp>
#include "energy.hpp"

namespace lur {

    class LurieObj : public Objective <double> {
    public:

        /**
         * Constructor
         * @param mm material model
         * @param potent potential
         */
        LurieObj(const lur::Energy& energy) : mEnergy(energy) {
            int n = energy.getModel().mNumLayers * 3;
            Objective<double>::setDim(n);
        }

        /**
         * Computes constraint value 
         * @param x parametes
         * @return constraint value
         */
        double func(const double* x) {
            double v = mEnergy.energy(x);            
            return v;
        }

       
    private:
        
        const lur::Energy& mEnergy;
    };

}

#endif	/* LUROBJ_HPP */

