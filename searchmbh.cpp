/**
 * Search for minimal energy with monotonic basin hopping methods
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
#include <problems/nlp/mbh/mbhboxcon.hpp>
#include <problems/nlp/mbh/staticpert.hpp>
#include <problems/nlp/mbh/adaptpert.hpp>
#include <util/box/boxutils.hpp>
#include <util/common/fileutils.hpp>
#include <util/spacefill/rndfill.hpp>
#include "lurobj.hpp"
#include "parsejson.hpp"
#include "genjson.hpp"
#include "ppenergy.hpp"
#include "tersoffparams.hpp"
#include "carbontersoff.hpp"
#include "tsofenergy.hpp"
#include "pairpotentials.hpp"


/**
 * Maximal local hops for basin Hopping
 */
const int H = 8;


/**
 * Vicinity size
 */
const double V = 2.0;

void initVicinity(int N, Box<double>& box) {
    for (int i = 0; i < N; i++) {
        box.mA[i] = -V;
        box.mB[i] = V;
    }
}

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
        //        std::cout << "n = " << n << "f = " << fval << ", gmin = " << gmin << ", xdiff = " << xdiff << ", fdiff = " << fdiff << "\n";
        if (n > 1000) {
            return true;
        } else
            return false;
    }
};

double x[100];

lur::MatModel mm;

double bv;

void output(int sig) {
    std::cout << "Found v = " << bv << "\n";
    VecUtils::vecPrint(3 * mm.mNumLayers, x);
    exit(-1);
}

/*
 * Testing the perturbers
 */
int main(int argc, char** argv) {

    if (argc != 2)
        BNB_ERROR_REPORT("Usage: searchmbh.exe json_file\n");
    std::string jsons;
    FileUtils::getStringFromFile(argv[1], jsons);
    lur::ParseJson::parseModelData(jsons, mm);
    double ev;
    lur::ParseJson::parseLatticeData(jsons, mm, ev, x);

#if 0
    lur::PotentialCutter pc(6, 0.5, lur::ljpotent);
    lur::PairPotentialEnergy enrg(mm, pc);
#endif
#if 0    
    // Morse
    lur::PairPotentialEnergy enrg(mm, morsepotent);
#endif
#if 1
    // Tersoff
    lur::TersoffParams tparam;
    lur::fillCarbonParametersTersoffOriginal(tparam);
    lur::TersoffUtils tutils(tparam);
    lur::TersoffEnergy enrg(mm, tutils);
    //enrg.setFixedAtoms(true);
#endif    

    const int N = mm.mNumLayers * 3;
    lur::LurieObj obj(enrg, mm);
    NumGradObjective<double> nobj(obj);
    nobj.setH(1E-8);

    Box<double> box(N);
    lur::ParseJson::parseBoxData(jsons, box);
    NlpProblem<double> prob;
    prob.mBox = box;
    prob.mObj = &obj;

    BBStopper stp;
    BBBoxDescent<double> locs(box, &stp);
    locs.getOptions().mHInit = 4;
    locs.getOptions().mDec = 0.5;
    locs.getOptions().mHLB = 1e-6;
    locs.getOptions().mInc = 1.75;
    locs.setObjective(&nobj);
    RndFill<double> rfill(box);

    Box<double> vicinity(N);
    initVicinity(N, vicinity);
#if 0    
    StaticPerturber<double> perturber(prob, vicinity);
#else
    AdaptPerturber<double> perturber(prob, vicinity, AdaptPerturber<double>::Params({1, .01, 2, 0.1, 1.1}));
#endif    
#if 0
    MBHBoxCon<double> mbhbc(prob, perturber, H);
#else    
    MBHBoxCon<double> mbhbc(prob, perturber, H, &locs);
#endif

    signal(SIGINT, output);

    bv = mbhbc.search(x);

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

