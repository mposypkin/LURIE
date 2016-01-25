
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
#include <problems/optlib/numgradobjective.hpp>
#include <problems/nlp/mbh/mbhboxcon.hpp>
#include <problems/nlp/mbh/staticpert.hpp>
#include <problems/nlp/mbh/adaptpert.hpp>
#include <problems/boxcon/boxconbnb/boxconbnb.hpp>
#include <problems/boxcon/boxconbnb/lipbounder.hpp>
#include <util/box/boxutils.hpp>
#include <problems/optlib/bbboxdesc.hpp>
#include <util/common/fileutils.hpp>
#include "lurobj.hpp"
#include "parsejson.hpp"
#include "genjson.hpp"
#include "ppenergy.hpp"
#include "tersoffparams.hpp"
#include "carbontersoff.hpp"
#include "tsofenergy.hpp"
#include "pairpotentials.hpp"

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
        if (n > 1000) {
            return true;
        } else
            return false;
    }
};

/*
 * Testing the perturbers
 */
int main(int argc, char** argv) {

    if (argc != 4)
        BNB_ERROR_REPORT("Usage: searchall.exe json_file log_for_incumbents log_for_allvalues\n");
    std::string jsons;
    FileUtils::getStringFromFile(argv[1], jsons);
    lur::MatModel mm;
    lur::ParseJson::parseModelData(jsons, mm);
    double x[mm.mNumLayers * 3];
    double y[mm.mNumLayers * 3];
    double ev;
    lur::ParseJson::parseLatticeData(jsons, mm, ev, x);

#if 0    
    // Lennard Jones
    lur::PairPotentialEnergy enrg(mm, ljpotent);
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

    Box<double> box(N);
    lur::ParseJson::parseBoxData(jsons, box);

    BBStopper stp;
    BBBoxDescent<double> locs(box, &stp);
    locs.getOptions().mHInit = 4;
    locs.getOptions().mDec = 0.5;
    locs.getOptions().mHLB = 1e-6;
    locs.getOptions().mInc = 1.75;
    locs.setObjective(&obj);


    BoxconProblem<double> bp;
    bp.mBox = box;
    bp.mObj = &obj;
    LipBounder<double> lb(bp.mObj);
    lb.setLipConst(96);
    BoxconBNB<double> bnb;
    bnb.addBounder(&lb);
    bnb.init(bp);

    double bestv = 0;
    int cnt = 0;
    time_t tstart = time(NULL);
    auto hndl = [&](BoxconBNB<double>::Solution& incum, BoxconBNB<double>::Sub & sub) {
        BoxUtils::getCenter(sub.mBox, y);
        double u = obj.func(y);
        //std::cout << "u = " << u << "\n";        
        // if ((u < 0.9 * incum.mValue) && ((cnt++ % 100) == 0)) {
        //if ((cnt++ % 100) == 0) {
        if (true) {
            std::cout << "Search from " << u << "\n";
            std::cout << "In box " << BoxUtils::toString(sub.mBox) << "\n";
            std::cout << " with radius " << BoxUtils::radius(sub.mBox) << "\n";

            double v;
            locs.search(y, &v);

            auto logval = [&] (const char* fname, double v) {
                time_t tcur = time(NULL) - tstart;
                std::ostringstream os;
                os << tcur << " " << v << "\n";
                FileUtils::updateFileWithContent(fname, os.str().c_str());
            };
            
            std::cout << "Found v = " << v << "\n";
            logval(argv[3], v);
            
            if (v < bestv) {
                bestv = v;
                std::cout << "Improved incumbent " << bestv << " from " << u << "\n";
                incum.mValue = v;
                VecUtils::vecPrint(N, y);
                VecUtils::vecCopy(N, y, x);
                logval(argv[2], v);
            }
        }
    };

    bnb.setHandler(hndl);
    long long steps = 100000;
    bnb.solve(steps);
#if 0
    std::cout << "The solution is " << bnb.getIncumbent().mValue << " found in " << steps << " steps \n";
    VecUtils::vecPrint(N, (double*) bnb.getIncumbent().mX);
    VecUtils::vecCopy(N, (double*) bnb.getIncumbent().mX, x);
#endif    



#if 0
    std::cout << "Searching in box " << BoxUtils::toString(box) << "\n";
    double v = mbhbc.search(x);
#endif

    //enrg.setFixedAtoms(true);
    double v;
    locs.search(x, &v);


    std::cout << "Found v = " << v << "\n";
    VecUtils::vecPrint(N, x);
    std::string json;
    json += "{\n";
    lur::GenJSON::genModel(mm, json);
    json += ", \n";
    lur::GenJSON::genBox(box, json);
    json += ", \n";
    lur::GenJSON::genLattice(mm, v, x, json);
    json += "\n}\n";
    std::cout << json << "\n";
    return 0;
}

