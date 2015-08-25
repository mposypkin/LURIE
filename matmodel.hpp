/* 
 * File:   matmodel.hpp
 * Author: medved
 * 
 * Layered material model definition
 * 
 * Created on August 21, 2015, 1:30 PM
 */

#ifndef MATMODEL_HPP
#define	MATMODEL_HPP



namespace lur {

    static const unsigned int maxLayers = 100;
    
    /**
     * Material model definition
     */
    struct MatModel {
        /**
         * Number of main layers
         */
        int mNumLayers;
        
        /**
         * The maximal radius of potential interaction
         */
        double mRadius;
               
        /**
         * The length of the material piece to model
         */
        double mLength;
        
        /**
         * Layer's atoms, each atom is identified by an integral number
         */
        int mLayersAtoms[maxLayers];
    };
}

#endif	/* MATMODEL_HPP */

