#define DRIVER_VERSION "0"

#include<iostream>
#include<fstream>
#include<sstream>
#include<unordered_map>
#include<unordered_set>
#include<cstdio>
#include<cstdlib>
#include<cstring>

#define ROOTING
#define NAME_MAPPING

typedef double score_t;

#include "landmark.hpp"
#include "algorithms.hpp"

using namespace std;

MetaAlgorithm meta;
TripartitionInitializer& tripInit = meta.tripInit;

vector<string>& names = meta.names;
unordered_map<string, int>& name2id = meta.name2id;

void addName(const string& name) {
    if (name2id.count(name) == 0) {
        name2id[name] = names.size();
        names.push_back(name);
    }
}

int main(int argc, char** argv){
    ARG.setProgramName("monster", "Morphological Oriented Novel Species Tree EstimatoR");
    
    meta.initialize(argc, argv);
    if (ARG.getStringArg("root") != "") addName(ARG.getStringArg("root"));

    ifstream fin(ARG.getStringArg("input"));
    string line;
    unordered_map<string, int> axisName2axis;
    int m = 0;
    vector<unordered_map<string, vector<double> > > axisInd2Array(3);
    while (getline(fin, line)){
        stringstream strin(line);
        double v;
        string indName, axisName;
        vector<double> arr;
        strin >> indName >> axisName;
        while (strin >> v) arr.push_back(v);
        if (m == 0) m = arr.size();
        if (axisName2axis.count(axisName) == 0){
            axisName2axis[axisName] = axisName2axis.size();
        }
        addName(meta.mappedname(indName));
        axisInd2Array[axisName2axis[axisName]][indName] = arr;
    }
    tripInit.M.resize(names.size(), vector<Vector>(m));
    for (auto &e: axisInd2Array[0]){
        string indName = e.first;
        int speciesID = name2id[meta.mappedname(indName)];
        for (int i = 0; i < m; i++){
            int j = (i == 0) ? m - 1 : i - 1;
            if (axisName2axis.size() == 2){
                Vector vd(axisInd2Array[0][indName][i], axisInd2Array[1][indName][i], 0, true);
                /*
                Vector vi(axisInd2Array[0][indName][i], axisInd2Array[1][indName][i], 0);
                Vector vj(axisInd2Array[0][indName][j], axisInd2Array[1][indName][j], 0);
                Vector vd = vi - vj;
                vd.normalize();
                */
                tripInit.M[speciesID][i] += vd;
            }
            if (axisName2axis.size() == 3){
                Vector vd(axisInd2Array[0][indName][i], axisInd2Array[1][indName][i], axisInd2Array[2][indName][i], true);
                /*
                Vector vi(axisInd2Array[0][indName][i], axisInd2Array[1][indName][i], axisInd2Array[2][indName][i]);
                Vector vj(axisInd2Array[0][indName][j], axisInd2Array[1][indName][j], axisInd2Array[2][indName][j]);
                Vector vd = vi - vj;
                vd.normalize();
                */
                tripInit.M[speciesID][i] += vd;
            }
        }
    }

    LOG << "#Dimensions: " << axisName2axis.size() << endl;
    LOG << "#Lanmarks: " << m << endl;
    LOG << "#Samples: " << axisInd2Array[0].size() << endl;
    
    auto res = meta.run();
    LOG << "Score: " << (double) res.first << endl;
    
    return 0;
}
