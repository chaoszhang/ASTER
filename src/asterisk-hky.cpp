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

void addName(const string &name){
    if (name2id.count(name) == 0){
        name2id[name] = names.size();
        names.push_back(name);
    }
}

void addBin(const vector<int> &geneTaxa, const vector<vector<char> > &seqs,
        const vector<array<score_t, 4> > &freq, const vector<array<int, 4> > &cnt){
    for (int j = 0; j < freq.size(); j++){
        if (cnt[j][0] + cnt[j][2] < 2 || cnt[j][1] + cnt[j][3] < 2) continue;
        for (int i = 0; i < seqs.size(); i++){
            int id = geneTaxa[i];
            tripInit.seq[id].push_back(seqs[i][j]);
        }
        tripInit.pi.push_back(freq[j]);
        tripInit.npos++;
    }
    for (int i = 0; i < tripInit.seq.size(); i++){
        while (tripInit.seq[i].size() < tripInit.npos) tripInit.seq[i].push_back(-1); 
    }
}


void formatGene(const vector<int> geneTaxa, const vector<string> geneSeqs){
	vector<vector<char> > seqs(geneSeqs.size());
    vector<array<int, 4> > cnt(geneSeqs[0].size());
    vector<array<score_t, 4> > freq(geneSeqs[0].size());
    int sum[4] = {}; 
    for (int i = 0; i < geneSeqs.size(); i++){
        for (int j = 0; j < geneSeqs[i].size(); j++){
            int x = -1;
            if (geneSeqs[i][j] == 'A' || geneSeqs[i][j] == 'a') x = 0;
            else if (geneSeqs[i][j] == 'C' || geneSeqs[i][j] == 'c') x = 1;
            else if (geneSeqs[i][j] == 'G' || geneSeqs[i][j] == 'g') x = 2;
            else if (geneSeqs[i][j] == 'T' || geneSeqs[i][j] == 't') x = 3;
            else continue;
            seqs[i].push_back(x);
            cnt[j][x]++;
            sum[x]++;
        }
    }
    int sum4 = sum[0] + sum[1] + sum[2] + sum[3];
    for (int j = 0; j < geneSeqs[0].size(); j++){
        score_t cnt4 = sum4 - (cnt[j][0] + cnt[j][1] + cnt[j][2] + cnt[j][3]);
        for (int x = 0; x < 4; x++) freq[j][x] = (sum[x] - cnt[j][x]) / cnt4;
    }

    addBin(geneTaxa, seqs, freq, cnt);
}

void readPhilip(istream &fin){
	int M, L;
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
		formatGene(geneTaxa, geneSeqs);
	}
}

void readFasta(string file){
	ifstream fin(file);
    vector<int> geneTaxa;
	vector<string> geneSeqs;
	string line;
	while (getline(fin, line)){
		if (line[0] == '>'){
			geneTaxa.push_back(name2id[formatName(line)]);
			geneSeqs.emplace_back();
		}
		else geneSeqs.back() += formatName(line);
	}
	formatGene(geneTaxa, geneSeqs);
}

string HELP_TEXT = R"V0G0N(-y  take one input in PHYLIP format instead of a list of inputs in FASTA format 
inputList: the path to a file containing a list of paths to input aligned gene files, one file per line
Gene files must be in FASTA format. The header line should be ">Species_Name".
)V0G0N";

int main(int argc, char** argv){
	ARG.setProgramName("asterisk-hky", "Accurate Species Tree EstimatoR from Independent Site Kernals");
	ARG.addStringArg('f', "format", "fasta", "Input file type, fasta: a txt file containing a list of FASTA gene files; phylip: a phylip file with all genes", true);
    ARG.addStringArg('d', "diskcover", "", "Full binary reference spcies tree in Newick format for disk covering method", true);
    ARG.addIntArg('b', "diskcoverbins", 10, "The number weight bins for disk covering method", true);
	ARG.addFlag('F', "formatphylip", "Indicating input file as a PHYLIP file containing all genes (`-f phylip`)", [&](){
		ARG.getStringArg("format") = "phylip";
	}, true);
	string mappingFile;
	meta.initialize(argc, argv, " -y", HELP_TEXT);
	tripInit.nThreads = meta.nThreads;

	bool phylip = ARG.getStringArg("format").compare("fasta");
	
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
		}
	}
    
    tripInit.seq.resize(names.size());
    tripInit.speciesMapper.resize(names.size());
    if (ARG.getStringArg("diskcover").compare("")){
        string s;
        ifstream fin(ARG.getStringArg("diskcover"));
        getline(fin, s);
        Tree dcTree(s, name2id);
        mt19937 eng;
        vector<int> nbps = {1};
        int num = 0, nBins = ARG.getIntArg("diskcoverbins");
        double logmin = log(1.0/12), logmax = log(1.0/(log2(1.0*names.size())*8-6));
        for (int i = 0; i < nBins - 1; i++){
            nbps.push_back(1 + names.size() * exp(logmin + (logmax - logmin) * i / (nBins - 2)));
        }
        for (int nbp: nbps){
            vector<vector<int> > bins = dcTree.diskCoveringNBin(nbp, eng);
            long long nq = 0;
            for (const vector<int> &bin: bins){
                if (bin.size() < 4) continue;
                nq += SeqUtils::choose4(bin.size());
                for (int i: bin) tripInit.speciesMapper[i].push_back(num);
                num++;
            }
            for (const vector<int> &bin: bins){
                tripInit.w.push_back(1.0 / nq);
            }
        }
    }
    else{
        tripInit.w.push_back(1.0 / SeqUtils::choose4(names.size()));
        for (int i = 0; i < names.size(); i++) tripInit.speciesMapper[i].push_back(0);
    }

	if (phylip) {
        ifstream fin(ARG.getStringArg("input"));
		readPhilip(fin);
	}
	else {
        ifstream listIn(ARG.getStringArg("input"));
		for (string file; getline(listIn, file);){
			readFasta(file);
		}
	}
	
    //cerr << "rateModifier: " << rateModifier << endl;
	cerr << "#Informative Sites: " << tripInit.npos << endl;
	
	score_t score = meta.run().first;
	cerr << "Score: " << (double) score << endl;
	return 0;
}
