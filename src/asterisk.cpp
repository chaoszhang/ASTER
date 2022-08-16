#define DRIVER_VERSION "0"

#include<iostream>
#include<fstream>
#include<unordered_map>
#include<unordered_set>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<algorithm>
#include<random>
#include<thread>
#include<mutex>

#define LARGE_DATA
#ifdef LARGE_DATA
typedef long double score_t;
typedef long long count_t;
#else
typedef double score_t;
typedef int count_t;
#endif

#include "argparser.hpp"
#include "sequtils.hpp"
#include "sequence-pair.hpp"
#include "algorithms.hpp"

using namespace std;

MetaAlgorithm meta;
TripartitionInitializer &tripInit = meta.tripInit;

vector<string> &names = meta.names;
unordered_map<string, int> &name2id = meta.name2id;

double MIN_MUTATION_RATE = 0.01;
int NUM_MUTATION_RATE_BINS = 10;
double MUTATION_RATE_BIN_RATIO = 2;
double MAX_MUTATION_RATE_RATIO = 1.1;

bool useDC = false;

void addName(const string &name){
    if (name2id.count(name) == 0){
        name2id[name] = names.size();
        names.push_back(name);
    }
}

void addBin(const vector<int> &geneTaxa, const AlignmentHot &a, const AlignmentHot &b,
        int goodLen, int geneID, const vector<double> &rates){
    vector<int> cntA(a.nSites()), cntB(b.nSites());
    for (int i = 0; i < a.nTaxa(); i++){
        for (int j = 0; j < a.nSites(); j++){
            cntA[j] += a[i][j];
            cntB[j] += b[i][j];
        }
    }
    long long evenA = 0, evenB = 0, evenAB = 0, even2 = 0;
    long long oddA = 0, oddB = 0, oddAB = 0, odd2 = 0;
    for (int j = 0; j < a.nSites(); j += 2) {
        evenA += cntA[j];
        evenB += cntB[j];
        evenAB += cntA[j] * cntB[j];
        even2 += (cntA[j] + cntB[j]) * (cntA[j] + cntB[j]);
    }
    for (int j = 1; j < a.nSites(); j += 2) {
        oddA += cntA[j];
        oddB += cntB[j];
        oddAB += cntA[j] * cntB[j];
        odd2 += (cntA[j] + cntB[j]) * (cntA[j] + cntB[j]);
    }
    vector<score_t> freq(goodLen);
    for (int j = 0; j < goodLen; j += 2) {
        score_t AB = (evenA - cntA[j]) * (evenB - cntB[j]) - (evenAB - cntA[j] * cntB[j]);
        score_t total = (evenA + evenB - cntA[j] - cntB[j]) * (evenA + evenB - cntA[j] - cntB[j])
            - (even2 - (cntA[j] + cntB[j]) * (cntA[j] + cntB[j]));
        freq[j] = AB / total;
    }
    for (int j = 1; j < goodLen; j += 2) {
        score_t AB = (oddA - cntA[j]) * (oddB - cntB[j]) - (oddAB - cntA[j] * cntB[j]);
        score_t total = (oddA + oddB - cntA[j] - cntB[j]) * (oddA + oddB - cntA[j] - cntB[j])
            - (odd2 - (cntA[j] + cntB[j]) * (cntA[j] + cntB[j]));
        freq[j] = AB / total;
    }
    int binID = 0;
    double rateNextBin = MIN_MUTATION_RATE * MUTATION_RATE_BIN_RATIO;
    for (int j = 0; j < goodLen; j++){
        while (useDC && binID != NUM_MUTATION_RATE_BINS - 1 && rates[j] > rateNextBin){
            binID++;
            rateNextBin *= MUTATION_RATE_BIN_RATIO;
        }
        for (int i = 0; i < geneTaxa.size(); i++){
            int id = geneTaxa[i];
            tripInit.trees[geneID][binID].seqX[id].push_back(a[i][j]);
            tripInit.trees[geneID][binID].seqY[id].push_back(b[i][j]);
        }
        tripInit.trees[geneID][binID].piXpiY.push_back(freq[j]);
        tripInit.trees[geneID][binID].npos++;
    }
    /*
    for (int i = 0; i < tripInit.seqX.size(); i++){
        while (tripInit.seqX[i].size() < tripInit.npos) {
            tripInit.seqX[i].push_back(0);
            tripInit.seqY[i].push_back(0);
        }
    }
    */
}


void formatGene(const vector<int> geneTaxa, const vector<string> geneSeqs, 
        vector<pair<double, int> > &mutationRate, int geneID){
	vector<int> goodSites, badSites;
    vector<double> rates;
    if (mutationRate.size()){
        sort(mutationRate.begin(), mutationRate.end());
        for (int i = 0; i < mutationRate.size(); i++){
            if (mutationRate[i].first < MIN_MUTATION_RATE || i == mutationRate.size() - 1
                    || mutationRate[i].first * MAX_MUTATION_RATE_RATIO < mutationRate[i+1].first)
                badSites.push_back(mutationRate[i].second);
            else {
                rates.push_back(mutationRate[i].first);
                rates.push_back(mutationRate[i].first);
                goodSites.push_back(mutationRate[i].second);
                i++;
                goodSites.push_back(mutationRate[i].second);
            }
        }
    }
    else{
        for (int i = 0; i < geneSeqs[0].size(); i++){
            if (i == geneSeqs[0].size() - 1) badSites.push_back(i);
            else {
                rates.push_back(MIN_MUTATION_RATE);
                rates.push_back(MIN_MUTATION_RATE);
                goodSites.push_back(mutationRate[i].second);
                i++;
                goodSites.push_back(mutationRate[i].second);
            }
        }
    }
    vector<string> orderedSeqs(geneSeqs.size());
    for (int i = 0; i < orderedSeqs.size(); i++){
        for (int j: goodSites) orderedSeqs[i].push_back(geneSeqs[i][j]);
        for (int j: badSites) orderedSeqs[i].push_back(geneSeqs[i][j]);
    }
    MSA msa(orderedSeqs);
    AlignmentHot a[7], b[7];
    a[0] = msa.A(); b[0] = msa.C() + msa.G() + msa.T();
    a[1] = msa.C(); b[1] = msa.A() + msa.G() + msa.T();
    a[2] = msa.G(); b[2] = msa.A() + msa.C() + msa.T();
    a[3] = msa.T(); b[3] = msa.A() + msa.C() + msa.G();
    a[4] = msa.A() + msa.C(); b[4] = msa.G() + msa.T();
    a[5] = msa.A() + msa.G(); b[5] = msa.C() + msa.T();
    a[6] = msa.A() + msa.T(); b[6] = msa.C() + msa.G();
    for (int r = 0; r < 7; r++){
        addBin(geneTaxa, a[r], b[r], goodSites.size(), geneID, rates);
    }
}

void readPhilip(istream &fin, vector<vector<pair<double, int> > > &mutationRate){
	int M, L, geneID = 0;
	vector<thread> thrds;
	while (fin >> M){
		fin >> L;
        vector<int> geneTaxa;
		vector<string> geneSeqs;
		for (int i = 0; i < M; i++){
			string s, line;
			fin >> s >> line;
			geneTaxa.push_back(name2id[s]);
			geneSeqs.push_back(line);
		}
		formatGene(geneTaxa, geneSeqs, mutationRate[geneID], geneID);
        geneID++;
	}
}

void readFasta(string file, vector<pair<double, int> > &mutationRate, int geneID){
	ifstream fin(file);
    vector<int> geneTaxa;
	vector<string> geneSeqs;
	string line;
	while (getline(fin, line)){
		if (line[0] == '>'){
			geneTaxa.push_back(name2id[SeqUtils::fastaFormatName(line)]);
			geneSeqs.emplace_back();
		}
		else geneSeqs.back() += SeqUtils::fastaFormatName(line);
	}
	formatGene(geneTaxa, geneSeqs, mutationRate, geneID);
}

string HELP_TEXT = R"V0G0N(-y  take one input in PHYLIP format instead of a list of inputs in FASTA format 
inputList: the path to a file containing a list of paths to input aligned gene files, one file per line
Gene files must be in FASTA format. The header line should be ">Species_Name".
)V0G0N";

int main(int argc, char** argv){
	ARG.setProgramName("asterisk", "Accurate Species Tree EstimatoR from Independent Site Kernals");
	ARG.addStringArg('f', "format", "fasta", "Input file type, fasta: a txt file containing a list of FASTA gene files; phylip: a phylip file with all genes", true);
    ARG.addStringArg('m', "mutation", "", "Substitution rate file from Iqtree if assumming heterogeneous rates", true);
    ARG.addStringArg('d', "diskcover", "", "Full binary reference spcies tree in Newick format for disk covering method", true);
    ARG.addIntArg('b', "diskcoverbins", 3, "The number weight bins for disk covering method", true);
	ARG.addFlag('F', "formatphylip", "Indicating input file as a PHYLIP file containing all genes (`-f phylip`)", [&](){
		ARG.getStringArg("format") = "phylip";
	}, true);
	string mappingFile;
	meta.initialize(argc, argv, " -y", HELP_TEXT);
	tripInit.nThreads = meta.nThreads;

	bool phylip = ARG.getStringArg("format").compare("fasta");
    int k = 0;

    if (phylip) {
        int M, L;
        ifstream fin(ARG.getStringArg("input"));
        while (fin >> M){
            fin >> L;
            for (int i = 0; i < M; i++){
                string s, line;
                fin >> s >> line;
                addName(s);
            }
            k++;
        }
	}
	else {
		ifstream listIn(ARG.getStringArg("input"));
		for (string file; getline(listIn, file);){
            ifstream fin(file);
            string line;
            while (getline(fin, line)){
                if (line[0] == '>') addName(SeqUtils::fastaFormatName(line));
            }
            k++;
		}
	}
    
    if (ARG.getStringArg("diskcover").compare("")){
        string s;
        ifstream fin(ARG.getStringArg("diskcover"));
        mt19937 eng;
        int nFold = ARG.getIntArg("diskcoverbins");
        useDC = true;
        tripInit.trees.resize(k, vector<TripartitionInitializer::Tree>(NUM_MUTATION_RATE_BINS, names.size()));
        for (int i = 0; i < k; i++){
            getline(fin, s);
            Tree dcTree(s, name2id);
            double rate = MIN_MUTATION_RATE;
            for (int j = 0; j < NUM_MUTATION_RATE_BINS; j++){
                vector<vector<int> > bins = dcTree.diskCovering(eng, rate, nFold);
                for (const vector<int> &bin: bins){
                    if (bin.size() < 4) continue;
                    for (int s: bin) tripInit.trees[i][j].speciesMapper[s].push_back(tripInit.trees[i][j].ntree);
                    tripInit.trees[i][j].ntree++;
                }
                rate *= MUTATION_RATE_BIN_RATIO;
            }
        }
    }
    else{
        tripInit.trees.resize(k, vector<TripartitionInitializer::Tree>(1, names.size()));
        for (int j = 0; j < k; j++){
            tripInit.trees[j][0].ntree = 1;
            for (int i = 0; i < names.size(); i++) tripInit.trees[j][0].speciesMapper[i].push_back(0);
        }
    }

    vector<vector<pair<double, int> > > mutationRate(k);
    if (ARG.getStringArg("mutation") != "") mutationRate = SeqUtils::iqtreeRateParser(ARG.getStringArg("mutation"));

	if (phylip) {
        ifstream fin(ARG.getStringArg("input"));
		readPhilip(fin, mutationRate);
	}
	else {
        int i = 0;
        ifstream listIn(ARG.getStringArg("input"));
		for (string file; getline(listIn, file);){
            readFasta(file, mutationRate[i], i);
            i++;
		}
	}
	
    //cerr << "rateModifier: " << rateModifier << endl;
	cerr << "#Genes: " << tripInit.trees.size() << endl;
	
	score_t score = meta.run().first;
	cerr << "Score: " << (double) score << endl;
	return 0;
}
