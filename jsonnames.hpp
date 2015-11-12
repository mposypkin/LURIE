/* 
 * File:   jsonnames.hpp
 * Author: medved
 * 
 * 
 *
 * Created on September 30, 2015, 5:29 PM
 */

#ifndef JSONNAMES_HPP
#define	JSONNAMES_HPP

namespace lur {

    namespace JsonNames {
        /**
         * Model definition
         */
        static const char modelname[] = "model";

        /**
         * number of layers
         */
        static const char numlayers[] = "nlay";


        /**
         *  Model's potentials range
         */
        static const char range[] = "range";

        /**
         *  Length of the modeled piece
         */
        static const char length[] = "length";

        /** 
         * Height of the modeled piece
         */
        static const char height[] = "height";

        /** 
         * atom types
         */
        static const char atoms[] = "atoms";
        
        /**
         * Lattice vector data
         */
        static const char lattice[] = "lattice";
        
        /**
         * Energy value
         */
        static const char evalue[] = "v";
        
        /**
         * Lattice vector
         */
        static const char lvector[] = "x";
        
        /**
         * Search bounding box
         */
        static const char box[] = "box";
        
        /**
         * Lower bound for the search box
         */
        static const char boxa[] = "a";
        
        /**
         * Upper bound for the search box
         */
        static const char boxb[] = "b";

    }
}


#endif	/* JSONNAMES_HPP */

