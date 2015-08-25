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
         * @param x parameters, defining the layers level, displacements and interatomic distances
         * @return the energy
         */
        static double energy(const MatModel& model, const double* x) {
            double E = 0;
            for (int i = 0; i < model.mNumLayers; i++) {
                E += layerEnergy(i, model, x);
            }
            return E;
        }
    };

    private:

    static double layerEnergy(int i, const MatModel& model, const double* x) {
        double h = x[3 * i];
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

    static double atomEnergy(int i, int j, const MatModel& model, const double* x) {
        for (int k = 0; k < model.mNumLayers; k++) {

        }
    }

    static getLayerHeight(int i, const MatModel& model, const double* x) {
        double h = x[model.mNumLayers - 1];
        int m = floor((double) i / (double) k);
        int j = i - m * model.mNumLayers;
        double r = h * m + x[j];
        return r;
    }
}

#endif	/* ENERGY_HPP */

