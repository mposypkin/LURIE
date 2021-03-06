/**
 * Search for minimal energy with gradient descent
 */

/* 
 * File:   testmbhboxcon.cpp
 * Author: medved
 *
 * Created on September 18, 2015, 8:44 AM
 */

#include <signal.h>
#include <problems/optlib/gradboxdesc.hpp>
#include <problems/optlib/bbboxdesc.hpp>
#include <problems/optlib/numgradobjective.hpp>
#include <util/box/boxutils.hpp>
#include <util/common/fileutils.hpp>
#include "lurobj.hpp"
#include "parsejson.hpp"
#include "genjson.hpp"
#include "ppenergy.hpp"
#include "tersoffparams.hpp"
#include "carbontersoff.hpp"
#include "tsofenergy.hpp"
#include "pairpotentials.hpp"



/**
 * Vicinity size
 */
const double V = .5;

class GbdStopper : public GradBoxDescent<double>::Stopper {
public:

    /**
     * Constructor
     * @param box problem box
     */
    GbdStopper(const Box<double>& box) : mBox(box) {
    }

    bool stopnow(double xdiff, double fdiff, double* x, double* grad, double fval, int n) {
        //        std::cout << "n = " << n << "f = " << fval << ", gnorm = " << gnorm << ", xdiff = " << xdiff << ", fdiff = " << fdiff << "\n";
        double gnorm = 0;
        double const eps = 0.00000001;
        for (int i = 0; i < mBox.mDim; i++) {
            if (BNBABS(x[i] - mBox.mA[i]) <= eps) {
                if (grad[i] < 0)
                    gnorm += grad[i] * grad[i];
            } else if (BNBABS(x[i] - mBox.mB[i]) <= eps) {
                if (grad[i] > 0)
                    gnorm += grad[i] * grad[i];
            } else {
                gnorm += grad[i] * grad[i];
            }
        }
        std::cout << "n = " << n << "f = " << fval << ", gnorm = " << gnorm << ", xdiff = " << xdiff << ", fdiff = " << fdiff << "\n";

        if (gnorm < 0.01) {
            return true;
        } else if (n > 100000) {
            return true;
        } else
            return false;
    }

private:
    const Box<double>& mBox;
};

class BBStopper : public BBBoxDescent<double>::Stopper {
public:

    /**
     * Returns true when the search should stop
     * @param xdiff difference between old and new x
     * @param fdiff difference between old and new f value
     * @param gmin gradient minimal component
     * @param fval function value
     * @param n current step number 
     */
    bool stopnow(double xdiff, double fdiff, double gmin, double fval, int n) {
        std::cout << "n = " << n << " gmin = " << gmin << ", fval = " << fval << "\n";
        //if(gmin > -1e-6)
            //return true;
        if (n > 10000) {
            return true;
        } else
            return false;
    }
};

/*
 * Testing the perturbers
 */

double x[100];
double bv;
lur::MatModel mm;

void output(int sig) {
    std::cout << "Found v = " << bv << "\n";
    VecUtils::vecPrint(3 * mm.mNumLayers, x);
    exit(-1);
}

int main(int argc, char** argv) {

    if (argc != 2)
        BNB_ERROR_REPORT("Usage: searchgdsc.exe json_file\n");
    std::string jsons;
    FileUtils::getStringFromFile(argv[1], jsons);

    lur::ParseJson::parseModelData(jsons, mm);
    double ev;
    lur::ParseJson::parseLatticeData(jsons, mm, ev, x);

#if 0    
    // Lennard Jones

    lur::PairPotentialEnergy enrg(mm, lur::ljpotent);
#endif
#if 0    
    // Cutted Lennard Jones
    lur::PairPotentialEnergy enrg(mm, ljcutpotent);
#endif
#if 0    
    // Morse
    lur::PairPotentialEnergy enrg(mm, lur::morsepotent);
#endif

#if 0
    lur::PotentialCutter pc(6, 0.5, lur::ljpotent);
    lur::PairPotentialEnergy enrg(mm, pc);
#endif

#if 1
    // Tersoff
    lur::TersoffParams tparam;
    lur::fillCarbonParametersTersoffOriginal(tparam);
    lur::TersoffUtils tutils(tparam);
    lur::TersoffEnergy enrg(mm, tutils);
    //enrg.setFixedAtoms(true);
#endif    

#if 1    
    // Sets the fixed atoms attribute
    // enrg.setFixedAtoms(true);
#endif    

    const int N = mm.mNumLayers * 3;
    lur::LurieObj obj(enrg, mm);
    NumGradObjective<double> nobj(obj);
    nobj.setH(1E-4);
    Box<double> box(N);
    lur::ParseJson::parseBoxData(jsons, box);
    std::cout << "Searching in box " << BoxUtils::toString(box) << "\n";

#if 0    
    GbdStopper stp(box);
    GradBoxDescent<double> locs(box, &stp);
    locs.getOptions().mGInit = .1;
#endif
#if 1   
    BBStopper stp;
    BBBoxDescent<double> locs(box, &stp);
    locs.getOptions().mHInit = 4;
    locs.getOptions().mDec = 0.5;
    locs.getOptions().mHLB = 1e-6;
    locs.getOptions().mInc = 1.75;
    //locs.getOptions().mOnlyCoordinateDescent = true;
#endif

    // TMP DEBUG
#if 0    
    double xxx[12] = {0.866173, 0.982091, 0.982755, 0.950314, 2.08087, 1.84751, 0.4, 2.99746, 1.86217, 0.865676, 1.49372, 0.974428};
    std::cout << nobj.func(xxx) << "\n";
    double ggg[12];
    nobj.grad(xxx, ggg);
    VecUtils::vecPrint(12, ggg);
    exit(0);
#endif    
    // TMP


    signal(SIGINT, output);

    VecUtils::vecPrint(N, x);
    locs.setObjective(&nobj);
    locs.search(x, &bv);

    std::cout << "Found v = " << bv << "\n";
    VecUtils::vecPrint(N, x);

    std::string json;
    json += "{\n";
    lur::GenJSON::genModel(mm, json);
    json += ", \n";
    lur::GenJSON::genBox(box, json);
    json += ", \n";
    lur::GenJSON::genLattice(mm, bv, x, json);
    json += "\n}\n";
    std::cout << json << "\n";

    return 0;
}


