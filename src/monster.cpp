#define DRIVER_VERSION "0"

#define ROOTING

#include<iostream>
#include<fstream>
#include<sstream>
#include<unordered_map>
#include<unordered_set>
#include<cstdio>
#include<cstdlib>
#include<cstring>

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
    int m;
    while (getline(fin, line)){
        stringstream strin(line);
        double v;
        string indName, axisName;
        vector<double> arr;
        strin >> indName >> axisName;
        while (strin >> v) arr.push_back(v);
        
    }

    exit(0);
    //auto res = meta.run();
    //LOG << "Score: " << (double) res.first << endl;
    
    return 0;
}
