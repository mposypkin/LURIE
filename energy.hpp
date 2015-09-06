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
#include "matmodel.hpp"

namespace lur {

    class ComputeEnergy {
    public:

        /**
         * Compute the energy of the material piece
         * @param model reference to the model defining the piece
         * @param x parameters, defining the layers height, displacements and interatomic distances
         * @return the energy
         */
        static double energy(const MatModel& model, const double* x) {
            double E = 0;
            /* Atoms are numbered from 0 */
            for (int i = 0; i < model.mNumLayers; i++) {
                E += layerEnergy(i, model, x);
            }
            return E;
        }
    };

    private:

    /**
     * Computes interaction energy for a layer
     * @param i layer's number
     * @param model model data
     * @param x layer's data
     * @return energy value
     */
    static double layerEnergy(int i, const MatModel& model, const double* x) {
        double h = (i == 0) ? 0 : x[3 * i];
        double d = x[3 * i + 1];
        double w = x[3 * i + 2];
        int j = 0;
        double E = 0;
        while (d < model.mLength) {
            E += atomEnergy(i, j, model, x);
            d += w;
            j++;
        }
        return E;
    }

    /**
     * Computes interaction energy for an atom (interactions with all neighbours are accounted)
     * @param i atoms' layer
     * @param j atoms' number in a layer
     * @param model model data
     * @param x layer's data
     * @return energy value
     */
    static double atomEnergy(int i, int j, const MatModel& model, const double* x) {
        double y = getLayerHeight(i);
        int k = model.mNumLayers;
        int ii = i % k;
        double myx = x[3 * i + 1] + j * x[3 * i + 2];
        double v = 0;
        for (int l = i + 1;; l++) {            
            h = getLayerHeight(l);
            if (h - y > model.mRadius)
                break;
            int ll = l % k;
            double d = x[3 * ll + 1];
            double w = x[3 * ll + 2];
            int tl = ceil((myx - R - d)/w);
            int tu = floor((myx - R + d)/w);
            for(int t = tl; t < tu; t ++) {
                v + interEnergy(i, j, l, t);
            }
        }
        return v;
    }

  
    
    /**
     * Computes interaction energy of two atoms 
     * @param i1 first atom layer
     * @param j1 first atom number
     * @param i2 second atom layer
     * @param j2 second atom number
     * @return interaction energy
     */
    static double interEnergy(int i1, int j1, int i2, int j2) {
        /** TMP **/
        return 1;
    }
}

#endif	/* ENERGY_HPP */

