#define DRIVER_VERSION "0"

#include<iostream>
#include<fstream>
#include<sstream>
#include<unordered_map>
#include<unordered_set>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<cmath>

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
    ifstream fin(ARG.getStringArg("input"));
    int m = 0, nSample = 0;
    vector<vector<double> > &M = tripInit.M, &cnt = tripInit.cnt;

    fin >> nSample >> m;
    if (ARG.getStringArg("root") != "") {
        addName(ARG.getStringArg("root"));
        M.emplace_back(m);
        cnt.emplace_back(m);
    }

    
    for (int t = 0; t < nSample; t++){
        string indName;
        fin >> indName;
        addName(meta.mappedname(indName));
        int i = name2id[meta.mappedname(indName)];
        if (i == M.size()){
            M.emplace_back(m);
            cnt.emplace_back(m);
        }
        for (int j = 0; j < m; j++){
            string s;
            fin >> s;
            if (s == "NAN" || s == "NaN" || s == "NA" || s == "nan" || s == "na") continue;
            try{
                double v = stod(s);
                if (!isnan(v)){
                    M[i][j] += v;
                    cnt[i][j]++;
                }
            }
            catch(...){

            }
        }
        
    }

    for (int i = 0; i < 10; i++){
        cerr << names[i];
        for (int j = 0; j < 20; j++){
            cerr << "\t" << M[i][j] << ":" << cnt[i][j];
        }
        cerr << endl;
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
