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


/**
 * Maximal number of MC tries
 */
const int ntries = 100;

/**
 * Maximal local hops
 */
const int H = 100;



/**
 * Vicinity size
 */
const double V = .5;

/*
void initbox(Box<double>& box) {
    int n = box.mDim;
    for (int i = 0; i < n; i++) {
        if (i % 3 == 0) {
            box.mA[i] = 0;
            box.mB[i] = 3;
        } else if (i % 3 == 1) {
            box.mA[i] = 0;
            box.mB[i] = 3;
        } else {
            box.mA[i] = 0.1;
            box.mB[i] = 3;
        }
    }
}
 */

void initVicinity(int N, Box<double>& box) {
    for (int i = 0; i < N; i++) {
        box.mA[i] = -V;
        box.mB[i] = V;
    }
}

double sqrpotent(int a1, int a2, double q) {
    return (q - 1)*(q - 1) - 1;
}

double ljpotent(int a1, int a2, double q) {
    BNB_ASSERT(q >= 0);
    if (q == 0)
        q = 0.00001;
    double u = q * q * q;
    double v = u * u;
    double p = 1. / v - 2. / u;
    return p;
}

double morsepotent(int a1, int a2, double q) {
    BNB_ASSERT(q >= 0);
    double r = sqrt(q);
    const double rho = 10;
    double E = exp(rho * (1 - r));
    double p = E * (E - 2);
    return p;
}

double ljcutpotent(int a1, int a2, double q) {
    BNB_ASSERT(q >= 0);
    double p;
    if (q > 9) {
        p = 0;
    } else {
        if (q == 0)
            q = 0.00001;
        double u = q * q * q;
        double v = u * u;
        p = 1. / v - 2. / u;
        p -= ljpotent(a1, a2, 9);
    }
    return p;
}

class GbdStopper : public GradBoxDescent<double>::Stopper {
public:

    bool stopnow(double xdiff, double fdiff, double gnorm, double fval, int n) {
        //std::cout << "n = " << n << "f = " << fval << ", gnorm = " << gnorm << ", xdiff = " << xdiff << ", fdiff = " << fdiff << "\n";
        if (gnorm < 0.1) {
            return true;
        } else if (xdiff < 0.01) {
            return true;
        } else if (n > 16) {
            return true;
        } else
            return false;
    }
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
//        std::cout << "n = " << n << "f = " << fval << ", gmin = " << gmin << ", xdiff = " << xdiff << ", fdiff = " << fdiff << "\n";
        if (n > 100) {
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
    double y[mm.mNumLayers * 3];
    double ev;
    lur::ParseJson::parseLatticeData(jsons, mm, ev, x);

#if 0    
    // Lennard Jones
    lur::PairPotentialEnergy enrg(mm, ljpotent);
#endif
#if 0
    // Lennard Jones
    lur::PairPotentialEnergy enrg(mm, ljcutpotent);
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
    Box<double> vicinity(N);
    initVicinity(N, vicinity);

#if 0    
    GbdStopper stp;
    GradBoxDescent<double> locs(box, &stp);
    locs.getOptions().mGInit = .01;
#endif    
#if 1   
    BBStopper stp;
    BBBoxDescent<double> locs(box, &stp);
    locs.getOptions().mHInit = 4;
    locs.getOptions().mDec = 0.5;
    locs.getOptions().mHLB = 1e-6;
    locs.getOptions().mInc = 1.75;
#endif

    locs.setObjective(&nobj);
    RndFill<double> rfill(box);


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

    bv = 0;
    for (int i = 0; i < ntries; i++) {
        double v;
        rfill.getPoint(y);
        v = obj.func(y);
        std::cout << "Iteration " << i;
        std::cout << " starting from v = " << v << "\n";
        v = mbhbc.search(y);
        std::cout << "MBH found " << v << "\n";
        if (v < bv) {
            std::cout << "Record update: v = " << v << ",\n";
            VecUtils::vecPrint(N, y);
            bv = v;
            VecUtils::vecCopy(N, y, x);
        }
    }



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


