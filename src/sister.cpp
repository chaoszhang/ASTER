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

#include "numerical.hpp"
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
    ARG.setProgramName("sister", "Synteny/Shape Inspired Species Tree EstimatoR");
    
    meta.initialize(argc, argv);
    if (ARG.getStringArg("root") != "") addName(ARG.getStringArg("root"));

    ifstream fin(ARG.getStringArg("input"));
    string line;
    int m = 0, nSample = 0;
    vector<vector<double> > &M = tripInit.M, &cnt = tripInit.cnt;

    while (getline(fin, line)){
        nSample++;
        stringstream strin(line);
        double v;
        string indName;
        vector<double> arr;
        strin >> indName;
        while (strin >> v) arr.push_back(v);
        if (m == 0) m = arr.size();
        addName(meta.mappedname(indName));
        int i = name2id[meta.mappedname(indName)];
        if (i == M.size()){
            M.emplace_back(m);
            cnt.emplace_back(m);
        }
        for (int j = 0; j < m; j++){
            if (arr[j] < 0) continue;
            M[i][j] += arr[j];
            cnt[i][j]++;
        }
    }
    
    for (int j = 0; j < M[0].size(); j++){
        double sumM = 0, sumCnt = 0;
        for (int i = 0; i < M.size(); i++){
            sumM += M[i][j];
            sumCnt += cnt[i][j];
        }
        double avg = (sumCnt > 0) ? sumM / sumCnt : 0;
        for (int i = 0; i < M.size(); i++){
            M[i][j] -= avg * cnt[i][j];
        }
    }
    
    tripInit.nThreads = meta.nThreads;

    LOG << "#Measurements: " << m << endl;
    LOG << "#Samples: " << nSample << endl;
    
    auto res = meta.run();
    LOG << "Score: " << (double) res.first << endl;
    
    return 0;
}
