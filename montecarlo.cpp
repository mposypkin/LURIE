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
#include <limits>
#include <sstream>
#include <problems/optlib/gradboxdesc.hpp>
#include <problems/optlib/bbboxdesc.hpp>
#include <problems/optlib/numgradobjective.hpp>
#include <util/box/boxutils.hpp>
#include <util/common/fileutils.hpp>
#include <util/spacefill/rndfill.hpp>
#include <problems/nlp/mbh/mbhboxcon.hpp>
#include <problems/nlp/mbh/staticpert.hpp>
#include <problems/nlp/mbh/adaptpert.hpp>
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

        VecUtils::vecPrint(mBox.mDim, x);
        VecUtils::vecPrint(mBox.mDim, grad);

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
        //std::cout << "n = " << n << " gmin = " << gmin << "\n";
        //if(gmin > -1e-6)
         //   return true;
        if (n > 1000) {
            return true;
        } else
            return false;
    }
};

double x[100];
double bv;
lur::MatModel mm;

void output(int sig) {
    std::cout << "Found v = " << bv << "\n";
    VecUtils::vecPrint(3 * mm.mNumLayers, x);
    exit(-1);
}

int main(int argc, char** argv) {

    if (argc != 5)
        BNB_ERROR_REPORT("Usage: searchgdsc.exe json_file monte_carlo_tries log_for_incumbents log_for_values\n");
    std::string jsons;
    FileUtils::getStringFromFile(argv[1], jsons);
    lur::ParseJson::parseModelData(jsons, mm);
    double ev;
    lur::ParseJson::parseLatticeData(jsons, mm, ev, x);
    int ntries = atoi(argv[2]);

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

#if 0    
    // Sets the fixed atoms attribute
    enrg.setFixedAtoms(true);
#endif    

    const int N = mm.mNumLayers * 3;
    lur::LurieObj obj(enrg, mm);

    Box<double> box(N);
    lur::ParseJson::parseBoxData(jsons, box);
    NlpProblem<double> prob;
    prob.mBox = box;
    prob.mObj = &obj;
    std::cout << "Searching in box " << BoxUtils::toString(box) << "\n";

    BBStopper stp;
    BBBoxDescent<double> locs(box, &stp);
    locs.getOptions().mHInit = 4;
    locs.getOptions().mDec = 0.5;
    locs.getOptions().mHLB = 1e-6;
    locs.getOptions().mInc = 1.75;
    //locs.getOptions().mOnlyCoordinateDescent = true;
    locs.setObjective(&obj);



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

    double y[N];
    double cv = std::numeric_limits<double>::max();
    ;
    int cnt = 0;

    time_t tstart = time(NULL);

    auto process = [&] () {
        auto logval = [&] (const char* fname, double v) {
            time_t tcur = time(NULL) - tstart;
            std::ostringstream os;
            os << tcur << " " << v << "\n";
            FileUtils::updateFileWithContent(fname, os.str().c_str());
        };
        cnt++;
        logval(argv[4], cv);
        std::cout << "Monte Carlo Iteration = " << cnt << "\n";
        if (cv < bv) {
            logval(argv[3], cv);
            std::cout << "Improved v = " << cv << "\n";
            VecUtils::vecPrint(N, y);
            bv = cv;
            VecUtils::vecCopy(N, y, x);
        }
    };

    signal(SIGINT, output);


    for (int i = 0; i < ntries; i++) {
        rfill.getPoint(y);
#if 1        
        locs.search(y, &cv);
#endif
#if 0
        cv = mbhbc.search(y);
#endif        
        process();
    }

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



