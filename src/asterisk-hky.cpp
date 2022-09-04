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
#include "sequence.hpp"
#include "algorithms.hpp"

using namespace std;

MetaAlgorithm meta;
TripartitionInitializer &tripInit = meta.tripInit;

vector<string> &names = meta.names;
unordered_map<string, int> &name2id = meta.name2id;

double MIN_MUTATION_RATE = 0.01;
int NUM_MUTATION_RATE_BINS = 10;
double MUTATION_RATE_BIN_RATIO = 2;

void addName(const string &name){
    if (name2id.count(name) == 0){
        name2id[name] = names.size();
        names.push_back(name);
    }
}

void addBin(const vector<int> &geneTaxa, const vector<vector<char> > &seqs, const vector<array<float, 4> > &freq, 
        const vector<array<int, 4> > &cnt, int geneID, vector<pair<double, int> > &mutationRate){
    int binID = 0;
    double rateNextBin = MIN_MUTATION_RATE * MUTATION_RATE_BIN_RATIO;
    vector<pair<double, int> > rate;
    if (mutationRate.size()){
        for (const auto &e: mutationRate) rate.push_back(e);
        sort(rate.begin(), rate.end());
    }
    else {
        for (int j = 0; j < freq.size(); j++) rate.push_back({(double) MIN_MUTATION_RATE, (int) j});
    }
    
    for (const auto &e: rate){
        while (mutationRate.size() && binID != NUM_MUTATION_RATE_BINS - 1 && e.first > rateNextBin){ //* 0.125 > rateNextBin){
            binID++;
            rateNextBin *= MUTATION_RATE_BIN_RATIO;
        }
        int j = e.second;
        if (cnt[j][0] + cnt[j][2] < 2 || cnt[j][1] + cnt[j][3] < 2) continue;
        //cerr << e.first << " ";
        for (int i = 0; i < seqs.size(); i++){
            int id = geneTaxa[i];
            tripInit.trees[geneID][binID].seq[id].push_back(seqs[i][j]);
        }
        tripInit.trees[geneID][binID].pi.push_back(freq[j]);
        tripInit.trees[geneID][binID].npos++;
    }
    /*
    for (int i = 0; i < tripInit.seq.size(); i++){
        while (tripInit.seq[i].size() < tripInit.npos) tripInit.seq[i].push_back(-1); 
    }
    */
}


void formatGene(const vector<int> geneTaxa, const vector<string> geneSeqs, 
        vector<pair<double, int> > &mutationRate, int geneID){
	vector<vector<char> > seqs(geneSeqs.size());
    vector<array<int, 4> > cnt(geneSeqs[0].size());
    vector<array<float, 4> > freq(geneSeqs[0].size());
    int sum[4] = {}; 
    for (int i = 0; i < geneSeqs.size(); i++){
        for (int j = 0; j < geneSeqs[i].size(); j++){
            int x = -1;
            if (geneSeqs[i][j] == 'A' || geneSeqs[i][j] == 'a') x = 0;
            else if (geneSeqs[i][j] == 'C' || geneSeqs[i][j] == 'c') x = 1;
            else if (geneSeqs[i][j] == 'G' || geneSeqs[i][j] == 'g') x = 2;
            else if (geneSeqs[i][j] == 'T' || geneSeqs[i][j] == 't') x = 3;
            seqs[i].push_back(x);
            if (x == -1) continue;
            cnt[j][x]++;
            sum[x]++;
        }
    }
    int sum4 = sum[0] + sum[1] + sum[2] + sum[3];
    for (int j = 0; j < geneSeqs[0].size(); j++){
        score_t cnt4 = sum4 - (cnt[j][0] + cnt[j][1] + cnt[j][2] + cnt[j][3]);
        for (int x = 0; x < 4; x++) freq[j][x] = (sum[x] - cnt[j][x]) / cnt4;
    }

    addBin(geneTaxa, seqs, freq, cnt, geneID, mutationRate);
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
		else geneSeqs.back() += SeqUtils::fastaFormatRead(line);
	}
	formatGene(geneTaxa, geneSeqs, mutationRate, geneID);
}

string HELP_TEXT = R"V0G0N(-y  take one input in PHYLIP format instead of a list of inputs in FASTA format 
inputList: the path to a file containing a list of paths to input aligned gene files, one file per line
Gene files must be in FASTA format. The header line should be ">Species_Name".
)V0G0N";

int main(int argc, char** argv){
	ARG.setProgramName("asterisk-hky", "Accurate Species Tree EstimatoR from Independent Site Kernals");
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
	mt19937 eng;
    eng.seed(ARG.getIntArg("seed"));

    if (phylip) {
        {
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
        
        if (ARG.getStringArg("diskcover").compare("")){
            string s;
            ifstream fin(ARG.getStringArg("diskcover"));
            int nFold = ARG.getIntArg("diskcoverbins");
            //useDC = true;
            tripInit.trees.resize(k, vector<TripartitionInitializer::Tree>(NUM_MUTATION_RATE_BINS, names.size()));
            for (int i = 0; i < k; i++){
                getline(fin, s);
                Tree dcTree(s, name2id);
                double rate = MIN_MUTATION_RATE * sqrt(MUTATION_RATE_BIN_RATIO);
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
        
        {
            ifstream fin(ARG.getStringArg("input"));
            readPhilip(fin, mutationRate);
        }
	}
	else {
        {
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

        int t = ARG.getIntArg("thread");
        if (ARG.getStringArg("diskcover").compare("")){
            string s;
            ifstream listIn(ARG.getStringArg("diskcover"));
            int nFold = ARG.getIntArg("diskcoverbins");
            //useDC = true;
            tripInit.trees.resize(k, vector<TripartitionInitializer::Tree>(NUM_MUTATION_RATE_BINS, names.size()));
            for (int i = 0; i < k; i++){
                string file;
                getline(listIn, file);
                ifstream fin(file);
                getline(fin, s);
                Tree dcTree(s, name2id);
                double rate = MIN_MUTATION_RATE * sqrt(MUTATION_RATE_BIN_RATIO);
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
            tripInit.trees.resize(t, vector<TripartitionInitializer::Tree>(1, names.size()));
            for (int j = 0; j < t; j++){
                tripInit.trees[j][0].ntree = 1;
                for (int i = 0; i < names.size(); i++) tripInit.trees[j][0].speciesMapper[i].push_back(0);
            }
        }
        
        if (ARG.getStringArg("mutation") != "") {
            ifstream listIn(ARG.getStringArg("input"));
            ifstream rateIn(ARG.getStringArg("mutation"));
            for (int i = 0; i < k; i++){
                string file, rfile;
                getline(listIn, file);
                getline(rateIn, rfile);
                readFasta(file, SeqUtils::iqtreeRateParser(rfile)[0], i);
                cerr << file << endl;
            }
        }
        else {
            ifstream listIn(ARG.getStringArg("input"));
            int i = 0;
            vector<pair<double, int> > mutationRate;
            for (string file; getline(listIn, file);){
                readFasta(file, mutationRate, i % t);
                if (i % 100 == 0) cerr << file << endl;
                i++;
            }
        }
	}
	
    //cerr << "rateModifier: " << rateModifier << endl;
	cerr << "#Genes: " << tripInit.trees.size() << endl;
	
	score_t score = meta.run().first;
	cerr << "Score: " << (double) score << endl;
	return 0;
}
