
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
#include <util/common/fileutils.hpp>
#include "lurobj.hpp"
#include "parsejson.hpp"
#include "genjson.hpp"
#include "ppenergy.hpp"
#include "tersoffparams.hpp"
#include "carbontersoff.hpp"
#include "tsofenergy.hpp"


/**
 * Maximal local hops
 */
const int H = 1000;

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

class MyStopper : public GradBoxDescent<double>::Stopper {
public:

    bool stopnow(double xdiff, double fdiff, double gnorm, double fval, int n) {
        std::cout << "n = " << n << "f = " << fval << ", gnorm = " << gnorm << ", xdiff = " << xdiff << ", fdiff = " << fdiff << "\n";
        if (gnorm < 0.0001) {
            return true;
        } else if (xdiff < 0.000000001) {
            return true;
        } else if (n > 1000) {
            return true;
        } else
            return false;
    }
};

/*
 * Testing the perturbers
 */
int main(int argc, char** argv) {

    if (argc != 2)
        BNB_ERROR_REPORT("Usage: searchall.exe json_file\n");
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
    NumGradObjective<double> nobj(obj);
    nobj.setH(1E-8);

    Box<double> box(N);
    lur::ParseJson::parseBoxData(jsons, box);
    NlpProblem<double> prob;
    prob.mBox = box;
    prob.mObj = &obj;
    Box<double> vicinity(N);
    initVicinity(N, vicinity);

    MyStopper stp;
    GradBoxDescent<double> gbd(box, &stp);
    gbd.getOptions().mGInit = .01;
    gbd.setObjective(&nobj);

#if 0    
    StaticPerturber<double> perturber(prob, vicinity);
#else
    AdaptPerturber<double> perturber(prob, vicinity, AdaptPerturber<double>::Params({1, .01, 2, 0.1, 1.1}));
#endif    
#if 1
    MBHBoxCon<double> mbhbc(prob, perturber, H);
#else    
    MBHBoxCon<double> mbhbc(prob, perturber, H, &gbd);
#endif

    BoxconProblem<double> bp;
    bp.mBox = box;
    bp.mObj = &nobj;
    LipBounder<double> lb(bp.mObj);
    lb.setLipConst(64);
    BoxconBNB<double> bnb;
    bnb.addBounder(&lb);
    bnb.init(bp);

    double bestv = 0;
    int cnt = 0;
    auto hndl = [&](BoxconBNB<double>::Solution& incum, BoxconBNB<double>::Sub & sub) {
        BoxUtils::getCenter(sub.mBox, y);
        double u = obj.func(y);
        //std::cout << "u = " << u << "\n";        
        if ((u < 0.9 * incum.mValue) && ((cnt++ % 100) == 0)) {
            //std::cout << "Search from " << u << "\n";
            double v;
#if 0
            gbd.search(y, &v);
#endif            
#if 1
            v = mbhbc.search(y);
#endif            
            //std::cout << "Found v = " << v << "\n";
            if (v < bestv) {
                bestv = v;
                std::cout << "Improved incumbent " << bestv << " from " << u << "\n";
                VecUtils::vecPrint(N, y);
                VecUtils::vecCopy(N, y, x);
            }
        }
    };

    //bnb.setHandler(hndl);
    long long steps = 100000;
    bnb.solve(steps);
#if 1
    std::cout << "The solution is " << bnb.getIncumbent().mValue << " found in " << steps << " steps \n";
    VecUtils::vecPrint(N, (double*) bnb.getIncumbent().mX);
    VecUtils::vecCopy(N, (double*) bnb.getIncumbent().mX, x);
#endif    



#if 0
    std::cout << "Searching in box " << BoxUtils::toString(box) << "\n";
    double v = mbhbc.search(x);
#endif
    
#if 1
    //enrg.setFixedAtoms(true);
    double v;
    gbd.search(x, &v);
#endif    
    
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

