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
#include "matmodel.hpp"

namespace lur {

    class LurieObj : public Objective <double> {
    public:

        /**
         * Constructor
         * @param mm material model
         * @param potent potential
         */
        LurieObj(lur::Energy& energy, const MatModel& model) : mEnergy(energy) {
            int n = model.mNumLayers * 3;
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
        
        lur::Energy& mEnergy;
                
    };

}

#endif	/* LUROBJ_HPP */

