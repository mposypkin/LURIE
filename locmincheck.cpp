
#include <problems/optlib/locminchecker.hpp>
#include <util/box/boxutils.hpp>
#include <util/common/fileutils.hpp>
#include "lurobj.hpp"
#include "ppenergy.hpp"
#include "parsejson.hpp"
#include "genjson.hpp"

double sqrpotent(int a1, int a2, double q) {
    return (q - 1)*(q - 1) - 1;
}

double ljpotent(int a1, int a2, double q) {
    BNB_ASSERT(q >= 0);
    if (q == 0)
        q = 0.00001;
    double u = q * q * q;
    double v = u * u;
    double p = 1. / v - 1. / u;
    return p;
}

double morsepotent(int a1, int a2, double q) {
    BNB_ASSERT(q >= 0);
    double r = sqrt(q);
    const double rho = .1;
    double E = exp(rho * (1 - r));
    double p = E * (E - 2);
    return p;
}

void makeBox(double* x, Box<double>& box, double eps) {
    int n = box.mDim;
    for (int i = 0; i < n; i++) {
        box.mA[i] = x[i] - eps;
        box.mB[i] = x[i] + eps;
    }
}

/*
 * Testing the perturbers
 */
int main(int argc, char** argv) {

    if (argc != 2)
        BNB_ERROR_REPORT("Usage: locmincheck.exe json_file\n");
    std::string jsons;
    FileUtils::getStringFromFile(argv[1], jsons);
    lur::MatModel mm;
    lur::ParseJson::parseModelData(jsons, mm);
    int n = mm.mNumLayers * 3;
    double x[n], y[n];
    double v;
    lur::ParseJson::parseLatticeData(jsons, mm, v, x);
    lur::PairPotentialEnergy enrg(mm, morsepotent);
    lur::LurieObj obj(enrg, mm);
    LocalMinChecker<double> lminchk(&obj);
    Box<double> box(n);
    makeBox(x, box, 0.001);
    int nt = 10000;
    bool fl = lminchk.check(x, y, nt, box);
    if (fl) {
        double vold = obj.func(x);
        std::cout << "Local minimum confirmed, v = " << vold << " \n";
    } else {
        std::cout << "Local minimum not confirmed\n";
        double vold = obj.func(x);
        double vnew = obj.func(y);
        std::cout << "Improved from " << vold << " to " << vnew << "\n";
    }

    return 0;
}
