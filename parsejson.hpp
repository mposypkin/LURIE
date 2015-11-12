/* 
 * File:   parsejson.hpp
 * Author: medved
 *
 * Parses JSON model definition
 * 
 * Created on September 30, 2015, 5:31 PM
 */

#ifndef PARSEJSON_HPP
#define	PARSEJSON_HPP

#include <string>
#include <libjson/libjson.h>
#include <util/common/bnberrcheck.hpp>
#include <util/box/box.hpp>
#include "matmodel.hpp"
#include "jsonnames.hpp"

namespace lur {

    class ParseJson {
    public:

        /**
         * Extracts JSON model from a string input
         * @param input JSON string
         * @param model model to fill in
         */
        static void parseModelData(const std::string& input, MatModel& model) {
            JSONNode nd = libjson::parse(input);
            bool modelset = false;
            for (auto i = nd.begin(); i != nd.end(); i++) {
                if (i->name() == JsonNames::modelname) {
                    processModel(*i, model);
                    modelset = true;
                    break;
                }
            }
            if (!modelset)
                BNB_ERROR_REPORT("Model data missing\n");
        }

        /**
         * Parses lattice definition vector from JSON string
         * @param input JSON string
         * @param model material model
         * @param v energy value to store
         * @param x vector to fill in
         */
        static void parseLatticeData(const std::string& input, const MatModel& model, double& v, double* x) {
            JSONNode nd = libjson::parse(input);
            bool latticeset = false;
            for (auto i = nd.begin(); i != nd.end(); i++) {
                if (i->name() == JsonNames::lattice) {
                    processLattice(*i, model.mNumLayers * 3, v, x);
                    latticeset = true;
                    break;
                }
            }
            if (!latticeset)
                BNB_ERROR_REPORT("Lattice data missing\n");

        }

        /**
         * Parses box data
         * @param input input string with JSON description
         * @param box to fill in
         */
        static void parseBoxData(const std::string& input, Box<double>& box) {
            JSONNode nd = libjson::parse(input);
            bool boxset = false;
            for (auto i = nd.begin(); i != nd.end(); i++) {
                if (i->name() == JsonNames::box) {
                    processBox(*i, box);
                    boxset = true;
                    break;
                }
            }
            if (!boxset)
                BNB_ERROR_REPORT("Lattice data missing\n");

        }


    private:

        static void processModel(const JSONNode & nd, MatModel& model) {
            bool nlayset = false;
            bool rangeset = false;
            bool lengthset = false;
            bool heightset = false;
            bool atomsset = false;
            for (auto i = nd.begin(); i != nd.end(); i++) {
                if (i->name() == JsonNames::numlayers) {
                    model.mNumLayers = i->as_int();
                    nlayset = true;
                } else if (i->name() == JsonNames::range) {
                    model.mRadius = i->as_float();
                    rangeset = true;
                } else if (i->name() == JsonNames::length) {
                    model.mLength = i->as_float();
                    lengthset = true;
                } else if (i->name() == JsonNames::height) {
                    model.mHeight = i->as_float();
                    heightset = true;
                } else if (i->name() == JsonNames::atoms) {
                    BNB_ASSERT(nlayset);
                    readIntVector(*i, model.mNumLayers, model.mLayersAtoms);
                    atomsset = true;
                } else {
                    BNB_ERROR_REPORT("Illegal name on parsing model data");
                }
            }
            BNB_ASSERT(nlayset && rangeset && lengthset && heightset && atomsset);
        }

        static void processBox(const JSONNode & nd, Box<double> & box) {
            int n = box.mDim;
            for (auto i = nd.begin(); i != nd.end(); i++) {
                if (i->name() == JsonNames::boxa) {
                    readDoubleVector(*i, n, (double*) (box.mA));
                } else if (i->name() == JsonNames::boxb) {
                    readDoubleVector(*i, n, (double*) (box.mB));
                } else {
                    BNB_ERROR_REPORT("Unknown name while processing box description");
                }
            }
        }

        static void processLattice(const JSONNode & nd, int n, double& v, double* x) {
            for (auto i = nd.begin(); i != nd.end(); i++) {
                if (i->name() == JsonNames::evalue) {
                    v = i->as_float();
                } else if (i->name() == JsonNames::lvector) {
                    readDoubleVector(*i, n, x);
                } else {
                    BNB_ERROR_REPORT("Unknown name while processing lattice description");
                }
            }

        }

        static void readIntVector(const JSONNode& nd, int vecsz, int * x) {
            int k = 0;
            for (auto i = nd.begin(); i != nd.end(); i++) {
                int u = i->as_int();
                BNB_ASSERT(k < vecsz);
                x[k++] = u;
            }
            BNB_ASSERT(k == vecsz);
        }

        static void readDoubleVector(const JSONNode& nd, int vecsz, double * x) {
            int k = 0;
            for (auto i = nd.begin(); i != nd.end(); i++) {
                double u = i->as_float();
                BNB_ASSERT(k < vecsz);
                x[k++] = u;
            }
            BNB_ASSERT(k == vecsz);
        }

    };
}


#endif	/* PARSEJSON_HPP */

