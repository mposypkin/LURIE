/**
 * Generates points in a rectangular region
 */

/* 
 * File: genpoints.cpp
 * Author: medved
 *
 * Generates points in a box. Stores results in files
 *  
 * Created on September 18, 2015, 8:44 AM
 */

#include <string>
#include <sstream>
#include <problems/nlp/mbh/mbhboxcon.hpp>
#include <problems/nlp/mbh/staticpert.hpp>
#include <util/box/boxutils.hpp>
#include <util/common/fileutils.hpp>
#include <util/spacefill/rndfill.hpp>
#include "lurobj.hpp"
#include "ppenergy.hpp"
#include "parsejson.hpp"
#include "genjson.hpp"

double ljpotent(int a1, int a2, double q) {
    BNB_ASSERT(q >= 0);
    if (q == 0)
        q = 0.00001;
    double u = q * q * q;
    double v = u * u;
    double p = 1. / v - 2. / u;
    return p;
}

int main(int argc, char** argv) {
    int ntries = 100;
    if (argc != 2)
        BNB_ERROR_REPORT("Usage: genpoints.exe json_file\n");
    std::string jsons;
    FileUtils::getStringFromFile(argv[1], jsons);
    lur::MatModel mm;
    lur::ParseJson::parseModelData(jsons, mm);
    int n = mm.mNumLayers * 3;
    lur::PairPotentialEnergy enrg(mm, ljpotent);
    const int N = mm.mNumLayers * 3;
    lur::LurieObj obj(enrg, mm);
    Box<double> box(N);
    double x[N];
    lur::ParseJson::parseBoxData(jsons, box);
    std::cout << "Searching in box " << BoxUtils::toString(box) << "\n";
    RndFill<double> rfill(box);
    double v;
    auto output = [ & ] (int I){
        std::string json;
        json += "{\n";
        lur::GenJSON::genModel(mm, json);
        json += ", \n";
        lur::GenJSON::genBox(box, json);
        json += ", \n";
        lur::GenJSON::genLattice(mm, v, x, json);
        json += "\n}\n";
        //std::cout << json << "\n";
        std::ostringstream os;
        os << I;        
        FileUtils::updateFileWithContent(os.str().c_str(), json.c_str());
    };
    for (int i = 0; i < ntries; i++) {
        rfill.getPoint(x);
        v = obj.func(x);
        VecUtils::vecPrint(n, x);
        std::cout << "Found v = " << v << "\n";
        output(i);
    }

    return 0;
}

