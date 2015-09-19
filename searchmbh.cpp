/**
 * Search for minimal energy with monotonic basin hopping methods
 */

/* 
 * File:   testmbhboxcon.cpp
 * Author: medved
 *
 * Created on September 18, 2015, 8:44 AM
 */

#include <problems/nlp/mbh/mbhboxcon.hpp>
#include <problems/nlp/mbh/staticpert.hpp>
#include <util/box/boxutils.hpp>
#include "lurobj.hpp"


/**
 * Maximal local hops
 */
const int H = 1000000;

/**
 * Vicinity size
 */
const double V = .5;

void initbox(Box<double>& box) {
    int n = box.mDim;
    for (int i = 0; i < n; i++) {
        if (i % 3 == 0) {
            box.mA[i] = 0;
            box.mB[i] = 2;
        } else if(i % 3 == 1) {
            box.mA[i] = 0;
            box.mB[i] = 2;            
        } else {
            box.mA[i] = 0.1;
            box.mB[i] = 2;            
        }
    }
}

void initVicinity(int N, Box<double>& box) {
    for (int i = 0; i < N; i++) {
        box.mA[i] = -V;
        box.mB[i] = V;
    }
}

void initx(int N, double* x) {
    for (int i = 0; i < N; i++) {
        x[i] = 1.5;
    }
}

double sqrpotent(int a1, int a2, double q) {
    return (q - 1)*(q - 1) - 1;
}


double ljpotent(int a1, int a2, double q) {
    BNB_ASSERT(q >= 0);
    if(q == 0)
        q = 0.00001;
    double u = q * q * q;
    double v = u * u;
    double p = 1./v - 1./u;
    return p;
}

/*
 * Testing the perturbers
 */
int main(int argc, char** argv) {
    lur::MatModel mm;
    mm.mNumLayers = 4;
    mm.mLength = 5.5;
    mm.mHeight = 6.5;
    mm.mRadius = 1.9;
//    lur::Energy enrg(mm, sqrpotent);
    lur::Energy enrg(mm, ljpotent);

    const int N = enrg.getModel().mNumLayers * 3;
    lur::LurieObj obj(enrg);
    Box<double> box(N);
    initbox(box);
    NlpProblem<double> prob;
    prob.mBox = box;
    prob.mObj = &obj;
    Box<double> vicinity(N);
    initVicinity(N, vicinity);
    StaticPerturber<double> perturber(prob, vicinity);
    MBHBoxCon<double> mbhbc(prob, perturber, H);
    double x[N];
    initx(N, x);

    std::cout << "Searching in box " << BoxUtils::toString(box) << "\n";
    double v = mbhbc.search(x);

    std::cout << "Found v = " << v << "\n";
    VecUtils::vecPrint(N, x);

    return 0;
}

