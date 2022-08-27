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

//#define LARGE_DATA
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
double MAX_MUTATION_RATE_RATIO = 1.02;
int MAX_SORT_BIN_SIZE = 1000;
int UNIFORM_RATE_STRIDE = 10;

bool useDC = false;
bool forceNeibourPairing = true;

void addName(const string &name){
    if (name2id.count(name) == 0){
        name2id[name] = names.size();
        names.push_back(name);
    }
}

void addBin(const vector<int> &geneTaxa, const AlignmentHot &a, const AlignmentHot &b,
        int goodLen, int geneID, const vector<double> &rates, double multiplier){
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
    for (int j = 0; j < goodLen; j += 2){
        while (useDC && binID != NUM_MUTATION_RATE_BINS - 1 && rates[j] * multiplier > rateNextBin){
            binID++;
            rateNextBin *= MUTATION_RATE_BIN_RATIO;
        }
        //cerr << cntA[j] << ' ' <<  cntB[j] << ' ' <<  cntA[j+1] << ' ' <<  cntB[j+1] << ' ' << endl;
        if ((cntA[j] <= 1 || cntB[j] <= 1) && (cntA[j + 1] == 0 || cntB[j + 1] == 0)) continue;
        if ((cntA[j] == 0 || cntB[j] == 0) && (cntA[j + 1] <= 1 || cntB[j + 1] <= 1)) continue;
        for (int i = 0; i < geneTaxa.size(); i++){
            int id = geneTaxa[i];
            tripInit.trees[geneID][binID].seqX[id].push_back(a[i][j]);
            tripInit.trees[geneID][binID].seqY[id].push_back(b[i][j]);
            tripInit.trees[geneID][binID].seqX[id].push_back(a[i][j + 1]);
            tripInit.trees[geneID][binID].seqY[id].push_back(b[i][j + 1]);
        }
        tripInit.trees[geneID][binID].piXpiY.push_back(freq[j]);
        tripInit.trees[geneID][binID].piXpiY.push_back(freq[j + 1]);
        tripInit.trees[geneID][binID].npos += 2;
    }
    /*
    for (int i = 0; i < tripInit.trees[geneID][binID].seqX.size(); i++){
        while (tripInit.trees[geneID][binID].seqX[i].size() < tripInit.trees[geneID][binID].npos) {
            tripInit.trees[geneID][binID].seqX[i].push_back(0);
            tripInit.trees[geneID][binID].seqY[i].push_back(0);
        }
    }
    */
}

void formatGenePartFourBin(const vector<int> &geneTaxa, const vector<string> &geneSeqs, 
        pair<double, int>* mutationRate, int mutationRateLen, int geneID, const vector<int> &alleleCnt){
	vector<int> goodSites, badSites;
    vector<double> rates;
    int usefulLen = (mutationRateLen / (UNIFORM_RATE_STRIDE * 4)) * (UNIFORM_RATE_STRIDE * 4);
    vector<tuple<double, int, int> > sorter;
    for (int i = 0; i < usefulLen; i += UNIFORM_RATE_STRIDE * 4){
        for (int j = i; j <= i + UNIFORM_RATE_STRIDE; j += UNIFORM_RATE_STRIDE){
            sort(&mutationRate[j + 0 * UNIFORM_RATE_STRIDE], &mutationRate[j + 1 * UNIFORM_RATE_STRIDE]);
            sort(&mutationRate[j + 2 * UNIFORM_RATE_STRIDE], &mutationRate[j + 3 * UNIFORM_RATE_STRIDE]);
            int a = j + 0 * UNIFORM_RATE_STRIDE, aEnd = j + 1 * UNIFORM_RATE_STRIDE;
            int b = j + 2 * UNIFORM_RATE_STRIDE, bEnd = j + 3 * UNIFORM_RATE_STRIDE;
            while (a < aEnd && b < bEnd){
                if (mutationRate[a].first * MAX_MUTATION_RATE_RATIO < mutationRate[b].first){ badSites.push_back(a); a++; continue; }
                if (mutationRate[b].first * MAX_MUTATION_RATE_RATIO < mutationRate[a].first){ badSites.push_back(b); b++; continue; }
                sorter.push_back({(double) mutationRate[a].first, (int) a, (int) b}); a++; b++;
            }
            while (a < aEnd) { badSites.push_back(a); a++; }
            while (b < bEnd) { badSites.push_back(b); b++; }
        }
    }
    for (int i = usefulLen; i < mutationRateLen; i++){
        badSites.push_back(i);
    }
    sort(sorter.begin(), sorter.end());
    for (const auto& e: sorter){
        rates.push_back(get<0>(e));
        rates.push_back(get<0>(e));
        goodSites.push_back(get<1>(e));
        goodSites.push_back(get<2>(e));
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
        double multiplier = (r == 4 || r == 6) ? 0.75 : 0.5;
        if (r == 4 || r == 6 || r == 5)
        addBin(geneTaxa, a[r], b[r], goodSites.size(), geneID, rates, multiplier);
    }
}

void formatGenePart(const vector<int> &geneTaxa, const vector<string> &geneSeqs, 
        pair<double, int>* mutationRate, int mutationRateLen, int geneID, const vector<int> &alleleCnt){
	vector<int> goodSites, badSites;
    vector<double> rates;
    formatGenePartFourBin(geneTaxa, geneSeqs, mutationRate, mutationRateLen, geneID, alleleCnt); return;
    if (mutationRateLen) {
        if (forceNeibourPairing){
            vector<tuple<double, int, int> > sorter;
            for (int i = 0; i < mutationRateLen; i += UNIFORM_RATE_STRIDE * 2){
                for (int j = i; j < mutationRateLen && j < i + UNIFORM_RATE_STRIDE; j++){
                    if (j + UNIFORM_RATE_STRIDE >= mutationRateLen) badSites.push_back(j);
                    else if (mutationRateLen && (mutationRate[j].first * MAX_MUTATION_RATE_RATIO < mutationRate[j + UNIFORM_RATE_STRIDE].first
                            || mutationRate[j + UNIFORM_RATE_STRIDE].first * MAX_MUTATION_RATE_RATIO < mutationRate[j].first)){
                        badSites.push_back(j);
                        badSites.push_back(j + UNIFORM_RATE_STRIDE);
                    }
                    else {
                        sorter.push_back({(double) mutationRate[j].first, (int) j, j + UNIFORM_RATE_STRIDE});
                    }
                }
            }
            sort(sorter.begin(), sorter.end());
            for (const auto& e: sorter){
                rates.push_back(get<0>(e));
                rates.push_back(get<0>(e));
                goodSites.push_back(get<1>(e));
                goodSites.push_back(get<2>(e));
            }
        }
        else {
            sort(mutationRate, &mutationRate[mutationRateLen]);
            for (int i = 0; i < mutationRateLen; i++){ 
                if (alleleCnt[mutationRate[i].second] <= 1) badSites.push_back(mutationRate[i].second); else
                if (i == mutationRateLen - 1 || mutationRate[i].first < MIN_MUTATION_RATE
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
    }
    else{
        for (int i = 0; i < geneSeqs[0].size(); i += UNIFORM_RATE_STRIDE * 2){
            for (int j = i; j < geneSeqs[0].size() && j < i + UNIFORM_RATE_STRIDE; j++){
                if (j + UNIFORM_RATE_STRIDE >= geneSeqs[0].size()) badSites.push_back(j);
                else {
                    rates.push_back(MIN_MUTATION_RATE);
                    rates.push_back(MIN_MUTATION_RATE);
                    goodSites.push_back(j);
                    goodSites.push_back(j + UNIFORM_RATE_STRIDE);
                }
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
        double multiplier = (r == 4 || r == 6) ? 0.75 : 0.5;
        //if (r == 4 || r == 6 || r == 5)
        addBin(geneTaxa, a[r], b[r], goodSites.size(), geneID, rates, multiplier);
    }
}

void formatGene(const vector<int> &geneTaxa, const vector<string> &geneSeqs, 
        vector<pair<double, int> > &mutationRate, int geneID){
	MSA msa(geneSeqs);
    vector<int> alleleCnt = msa.alleleCnt();
    vector<tuple<double, int, int> > sorter;
    if (mutationRate.size()){
        int nBin = (mutationRate.size() - 1) / MAX_SORT_BIN_SIZE + 1;
        for (int b = 0; b < nBin; b++){
            int start = b * mutationRate.size() / nBin, end = (b + 1) * mutationRate.size() / nBin;
            formatGenePart(geneTaxa, geneSeqs, &mutationRate[start], end - start, geneID, alleleCnt);
        }
    }
    else formatGenePart(geneTaxa, geneSeqs, &mutationRate[0], 0, geneID, alleleCnt);
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
            mt19937 eng;
            eng.seed(ARG.getIntArg("seed"));
            int nFold = ARG.getIntArg("diskcoverbins");
            useDC = true;
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

        int t = ARG.getIntArg("thread") * 10;
        if (ARG.getStringArg("diskcover").compare("")){
            string s;
            ifstream listIn(ARG.getStringArg("diskcover"));
            mt19937 eng;
            eng.seed(ARG.getIntArg("seed"));
            int nFold = ARG.getIntArg("diskcoverbins");
            useDC = true;
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
                readFasta(file, SeqUtils::iqtreeRateParser(rfile)[0], useDC ? i : i % t);
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
