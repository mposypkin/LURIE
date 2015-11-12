/* 
 * File:   genjson.hpp
 * Author: medved
 *
 * Generates JSON data for a model and a lattice
 * 
 * Created on October 1, 2015, 9:50 AM
 */

#ifndef GENJSON_HPP
#define	GENJSON_HPP

#include <iostream>
#include <sstream>

#include "jsonnames.hpp"
#include "matmodel.hpp"

namespace lur {

    class GenJSON {
    public:

        /**
         * Generates JSON description for a model
         * @param mm model
         * @param ss output string to append description
         */
        static void genModel(const MatModel& mm, std::string& ss) {
            std::ostringstream os;
            os << "\"" << JsonNames::modelname << "\" : {\n";
            os << "\"" << JsonNames::numlayers << "\" : " << mm.mNumLayers << ", \n";
            os << "\"" << JsonNames::range << "\" : " << mm.mRadius << ", \n";
            os << "\"" << JsonNames::length << "\" : " << mm.mLength << ", \n";
            os << "\"" << JsonNames::height << "\" : " << mm.mHeight << ", \n";
            os << "\"" << JsonNames::atoms << "\" : [";
            for(int i = 0; i < mm.mNumLayers; i ++) {
                os << mm.mLayersAtoms[i] << ((i == (mm.mNumLayers - 1)) ? "]\n" : ", ");
            }
            os << "}";
            ss += os.str();
        }

        /**
         * Generates JSON description of lattice parameters
         * @param mm mathematical model
         * @param v energy vector
         * @param x lattice vector
         * @param ss output string
         */
        static void genLattice(const MatModel& mm, double v, const double* x, std::string& ss) {
            std::ostringstream os;
            os << "\"" << JsonNames::lattice << "\" : {\n";
            os << "\"" << JsonNames::evalue << "\" : " << v << ",\n"; 
            os << "\"" << JsonNames::lvector << "\" : ["; 
            int n = mm.mNumLayers * 3;
            for(int i = 0; i < n; i ++) {
                os << x[i] << ((i == (n - 1)) ? "]" : ", ");
            }
            os << "\n}";
            ss += os.str();
        }
        
        /**
         * Generates JSON description for a box
         * @param box source box
         * @param ss output string
         */
        static void genBox(const Box<double>& box, std::string& ss) {
            int n = box.mDim;
            std::ostringstream os;
            os << "\"" << JsonNames::box << "\" : {\n";
            os << "\"" << JsonNames::boxa << "\" : ["; 
            for(int i = 0; i < n; i ++) {
                os << box.mA[i] << ((i == (n - 1)) ? "],\n" : ", ");
            }
            os << "\"" << JsonNames::boxb << "\" : ["; 
            for(int i = 0; i < n; i ++) {
                os << box.mB[i] << ((i == (n - 1)) ? "]\n" : ", ");
            }
            os << "}";
            ss += os.str();          
            
        }
    };
}

#endif	/* GENJSON_HPP */

