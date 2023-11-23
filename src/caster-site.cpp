#define DRIVER_VERSION "2"

/* CHANGE LOG
 * 2: Adding basic branch length functionality
 * 1: Modified the logic in parsing FASTA names 
 */

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

#define ROOTING
#define NAME_MAPPING
#define BLOCK_BOOTSTRAP

//#define LARGE_DATA
#ifdef LARGE_DATA
typedef long double score_t;
typedef long long count_t;
#else
typedef double score_t;
typedef int count_t;
#endif

#include "sequence.hpp"
#include "algorithms.hpp"
#include "sequtils.hpp"

using namespace std;

int GENE_ID = 0;

string getFastaName(const string &fasta){
    int i = 0;
    string res;
    while (i < fasta.size() && fasta[i] == ' ') i++;
    if (i == fasta.size() || fasta[i] != '>'){
        cerr << "Error in parsing line '" << fasta << "' in FASTA format.\n";
        exit(0);
    }
    i++;
    while (i < fasta.size() && fasta[i] == ' ') i++;
    while (i < fasta.size() && fasta[i] != ' ') {
        res += fasta[i];
        i++;
    }
    if (res.size() == 0){
        cerr << "Error in parsing line '" << fasta << "' in FASTA format.\n";
        exit(0);
    }
    return res;
}

pair<string, string> resolveDiploid(const string &seq){
	pair<string, string> res;
	for (const char c: seq){
		switch (c) {
			case 'A': case 'a': res.first.push_back('A'); res.second.push_back('A'); break;
			case 'C': case 'c': res.first.push_back('C'); res.second.push_back('C'); break;
			case 'G': case 'g': res.first.push_back('G'); res.second.push_back('G'); break;
			case 'T': case 't': res.first.push_back('T'); res.second.push_back('T'); break;
			case 'U': case 'u': res.first.push_back('T'); res.second.push_back('T'); break;
			case 'M': case 'm': res.first.push_back('A'); res.second.push_back('C'); break;
			case 'R': case 'r': res.first.push_back('A'); res.second.push_back('G'); break;
			case 'W': case 'w': res.first.push_back('A'); res.second.push_back('T'); break;
			case 'S': case 's': res.first.push_back('C'); res.second.push_back('G'); break;
			case 'Y': case 'y': res.first.push_back('C'); res.second.push_back('T'); break;
			case 'K': case 'k': res.first.push_back('G'); res.second.push_back('T'); break;
			default: res.first.push_back('N'); res.second.push_back('N');
		}
	}
	return res;
}

struct Workflow {
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

    static int infoCount(const array<int, 4> &cnt){
        int result = 0;
        for (int j = 0; j < 4; j++) {
            if (cnt[j] >= 2) result++;
        }
        return result;
    }

    static array<int, 4> add(const array<int, 4> &a, const array<int, 4> &b){
        array<int, 4> result;
        for (int j = 0; j < 4; j++) {
            result[j] = a[j] + b[j];
        }
        return result;
    }

    static int sum(const array<int, 4> &cnt){
        int result = 0;
        for (int j = 0; j < 4; j++) {
            result += cnt[j];
        }
        return result;
    }

    void ind2kernal(vector<int> &result, const unordered_map<int, int> &kernals, BinaryTree &hierarchy, int cur, int kernal = -1){
        if (kernals.count(cur)) kernal = kernals.at(cur);
        if (hierarchy.taxon(cur) != -1){
            result[hierarchy.taxon(cur)] = kernal;
        }
        else {
            int lc = hierarchy.left(cur), rc = hierarchy.right(cur);
            ind2kernal(result, kernals, hierarchy, lc, kernal);
            ind2kernal(result, kernals, hierarchy, rc, kernal);
        }
    }

    void buildGeneSeq(TripartitionInitializer::Gene::Initializer &gene, const vector<int> &ind2species, size_t pos, size_t nSite, size_t offset){
        int nInd = ind2species.size();
        array<size_t, 4> cnt = {};
        for (int iInd = 0; iInd < nInd; iInd++) {
            size_t pSeq = pos + iInd * offset;
            for (size_t iSite = 0; iSite < nSite; iSite++) {
                tripInit.seq.add(pSeq + iSite, cnt);
            }
            gene.species2ind[ind2species[iInd]].push_back(iInd);
            gene.ind2seq[iInd] = pSeq;
        }
        double total = cnt[0] + cnt[1] + cnt[2] + cnt[3];
        for (int j = 0; j < 4; j++) {
            gene.pi[j] = cnt[j] / total;
        }
    }

    void formatGene(const vector<int> &ind2species, size_t pos, size_t nSite, size_t offset) {
        int nInd = ind2species.size(), nSpecies = names.size(), nKernal = nSite, nRep = 0;
        TripartitionInitializer::Gene::Initializer gene(nInd, nSpecies, nSite, nKernal, nRep);
        buildGeneSeq(gene, ind2species, pos, nSite, offset);
        tripInit.genes.emplace_back(gene);
    }

    void readFasta(const string file) {
        cerr << "Processing " << file << "... \n";
        vector<bool> keep;
        {
            vector<array<unsigned short, 4> > freq;
            ifstream fin(file);
            string name;
            getline(fin, name);
            while (name != "") {
                addName(meta.mappedname(getFastaName(name)));
                name = "";
                string seq, line;
                while (getline(fin, line)) {
                    if (line[0] == '>') { name = line; break; }
                    seq += line;
                }
                if (freq.size() != seq.size()) {
                    if (freq.size() == 0) freq.resize(seq.size());
                    else { cerr << "File '" << file << "' is ill-formated.\n"; exit(0); }
                }
                for (size_t i = 0; i < seq.size(); i++) {
                    switch (seq[i]) {
                        case 'A': case 'a': freq[i][0]++; break;
                        case 'C': case 'c': freq[i][1]++; break;
                        case 'G': case 'g': freq[i][2]++; break;
                        case 'T': case 't': freq[i][3]++; break;
                        case 'U': case 'u': freq[i][3]++; break;
                    }
                }
            }
            keep.resize(freq.size());
            for (size_t i = 0; i < keep.size(); i++) {
                int cnt = 0, cnt1 = 0;
                for (int j = 0; j < 4; j++) {
                    if (freq[i][j] >= 2) cnt++;
                    if (freq[i][j] == 1) cnt1++;
                }
                keep[i] = (cnt >= 2 || cnt1 >= 2);
            }
        }
        {
            ifstream fin(file);
            string name;
            getline(fin, name);
            vector<int> ind2species;
            size_t pos = tripInit.seq.len(), len = count(keep.begin(), keep.end(), true);
            while (name != "") {
                ind2species.push_back(name2id[meta.mappedname(getFastaName(name))]);
                name = "";
                string line, seq;
                while (getline(fin, line)) {
                    if (line[0] == '>') { name = line; break; }
                    seq += line;
                }
                for (int i = 0; i < keep.size(); i++){
                    if (keep[i]) tripInit.seq.append(seq[i]);
                }
            }
            int nChunk = (len + ARG.getIntArg("chunk") - 1) / ARG.getIntArg("chunk");
            for (int i = 0; i < nChunk; i++) {
                size_t s = i * len / nChunk, t = (i + 1) * len / nChunk;
                formatGene(ind2species, pos + s, t - s, len);
            }
        }
    }
	
	void readFastaAmbiguity(const string file) {
        cerr << "Processing " << file << "... \n";
        vector<bool> keep;
        {
            vector<array<unsigned short, 4> > freq;
            ifstream fin(file);
            string name;
            getline(fin, name);
            while (name != "") {
                addName(meta.mappedname(getFastaName(name)));
                name = "";
                string seq, line;
                while (getline(fin, line)) {
                    if (line[0] == '>') { name = line; break; }
                    seq += line;
                }
                if (freq.size() != seq.size()) {
                    if (freq.size() == 0) freq.resize(seq.size());
                    else { cerr << "File '" << file << "' is ill-formated.\n"; exit(0); }
                }
				pair<string, string> seqs = resolveDiploid(seq);
                for (size_t i = 0; i < seq.size(); i++) {
                    switch (seqs.first[i]) {
                        case 'A': case 'a': freq[i][0]++; break;
                        case 'C': case 'c': freq[i][1]++; break;
                        case 'G': case 'g': freq[i][2]++; break;
                        case 'T': case 't': freq[i][3]++; break;
                        case 'U': case 'u': freq[i][3]++; break;
                    }
					switch (seqs.second[i]) {
                        case 'A': case 'a': freq[i][0]++; break;
                        case 'C': case 'c': freq[i][1]++; break;
                        case 'G': case 'g': freq[i][2]++; break;
                        case 'T': case 't': freq[i][3]++; break;
                        case 'U': case 'u': freq[i][3]++; break;
                    }
                }
            }
            keep.resize(freq.size());
            for (size_t i = 0; i < keep.size(); i++) {
                int cnt = 0, cnt1 = 0;
                for (int j = 0; j < 4; j++) {
                    if (freq[i][j] >= 2) cnt++;
                    if (freq[i][j] == 1) cnt1++;
                }
                keep[i] = (cnt >= 2 || cnt1 >= 2);
            }
        }
        {
            ifstream fin(file);
            string name;
            getline(fin, name);
            vector<int> ind2species;
            size_t pos = tripInit.seq.len(), len = count(keep.begin(), keep.end(), true);
            while (name != "") {
                ind2species.push_back(name2id[meta.mappedname(getFastaName(name))]);
                ind2species.push_back(name2id[meta.mappedname(getFastaName(name))]);
                name = "";
                string line, seq;
                while (getline(fin, line)) {
                    if (line[0] == '>') { name = line; break; }
                    seq += line;
                }
				pair<string, string> seqs = resolveDiploid(seq);
                for (int i = 0; i < keep.size(); i++){
                    if (keep[i]) tripInit.seq.append(seqs.first[i]);
                }
				for (int i = 0; i < keep.size(); i++){
                    if (keep[i]) tripInit.seq.append(seqs.second[i]);
                }
            }
            int nChunk = (len + ARG.getIntArg("chunk") - 1) / ARG.getIntArg("chunk");
            for (int i = 0; i < nChunk; i++) {
                size_t s = i * len / nChunk, t = (i + 1) * len / nChunk;
                formatGene(ind2species, pos + s, t - s, len);
            }
        }
    }
	
    bool readPhylip(ifstream &fin, ifstream &fin2){
        size_t nTaxa, nSites;
        if (!(fin >> nTaxa)) return false;
        fin >> nSites;
        vector<bool> keep(nSites);
        {
            vector<array<unsigned short, 4> > freq(nSites);
            for (int i = 0; i < nTaxa; i++){
                string name, seq;
                fin >> name >> seq;
                if (seq.size() != nSites) { cerr << "The input is ill-formated.\n"; exit(0); }
                addName(meta.mappedname(name));
                for (int j = 0; j < seq.size(); j++) {
                    switch (seq[j]) {
                        case 'A': case 'a': freq[j][0]++; break;
                        case 'C': case 'c': freq[j][1]++; break;
                        case 'G': case 'g': freq[j][2]++; break;
                        case 'T': case 't': freq[j][3]++; break;
						case 'U': case 'u': freq[j][3]++; break;
                    }
                }
            }
            for (int i = 0; i < keep.size(); i++) {
                int cnt = 0, cnt1 = 0;
                for (int j = 0; j < 4; j++) {
                    if (freq[i][j] >= 2) cnt++;
                    if (freq[i][j] == 1) cnt1++;
                }
                keep[i] = (cnt >= 2 || cnt1 >= 2);
            }
        }
        {
            fin2 >> nTaxa >> nSites;
            size_t pos = tripInit.seq.len(), len = count(keep.begin(), keep.end(), true);
            vector<int> ind2species;
            for (int i = 0; i < nTaxa; i++){
                string name, seq;
                fin2 >> name >> seq;
                ind2species.push_back(name2id[meta.mappedname(name)]);
                for (size_t j = 0; j < keep.size(); j++){
                    if (keep[j]) tripInit.seq.append(seq[j]);
                }
            }
            int nChunk = (len + ARG.getIntArg("chunk") - 1) / ARG.getIntArg("chunk");
            for (int i = 0; i < nChunk; i++) {
                size_t s = i * len / nChunk, t = (i + 1) * len / nChunk;
                formatGene(ind2species, pos + s, t - s, len);
            }
        }
        return true;
    }
	
	bool readPhylipAmbiguity(ifstream &fin, ifstream &fin2){
        size_t nTaxa, nSites;
        if (!(fin >> nTaxa)) return false;
        fin >> nSites;
        vector<bool> keep(nSites);
        {
            vector<array<unsigned short, 4> > freq(nSites);
            for (int i = 0; i < nTaxa; i++){
                string name, seq;
                fin >> name >> seq;
                if (seq.size() != nSites) { cerr << "The input is ill-formated.\n"; exit(0); }
                addName(meta.mappedname(name));
				pair<string, string> seqs = resolveDiploid(seq);
                for (int j = 0; j < seq.size(); j++) {
                    switch (seqs.first[j]) {
                        case 'A': case 'a': freq[j][0]++; break;
                        case 'C': case 'c': freq[j][1]++; break;
                        case 'G': case 'g': freq[j][2]++; break;
                        case 'T': case 't': freq[j][3]++; break;
						case 'U': case 'u': freq[j][3]++; break;
                    }
					switch (seqs.second[j]) {
                        case 'A': case 'a': freq[j][0]++; break;
                        case 'C': case 'c': freq[j][1]++; break;
                        case 'G': case 'g': freq[j][2]++; break;
                        case 'T': case 't': freq[j][3]++; break;
						case 'U': case 'u': freq[j][3]++; break;
                    }
                }
            }
            for (int i = 0; i < keep.size(); i++) {
                int cnt = 0, cnt1 = 0;
                for (int j = 0; j < 4; j++) {
                    if (freq[i][j] >= 2) cnt++;
                    if (freq[i][j] == 1) cnt1++;
                }
                keep[i] = (cnt >= 2 || cnt1 >= 2);
            }
        }
        {
            fin2 >> nTaxa >> nSites;
            size_t pos = tripInit.seq.len(), len = count(keep.begin(), keep.end(), true);
            vector<int> ind2species;
            for (int i = 0; i < nTaxa; i++){
                string name, seq;
                fin2 >> name >> seq;
                ind2species.push_back(name2id[meta.mappedname(name)]);
                ind2species.push_back(name2id[meta.mappedname(name)]);
				pair<string, string> seqs = resolveDiploid(seq);
                for (size_t j = 0; j < keep.size(); j++){
                    if (keep[j]) tripInit.seq.append(seqs.first[j]);
                }
				for (size_t j = 0; j < keep.size(); j++){
                    if (keep[j]) tripInit.seq.append(seqs.second[j]);
                }
            }
            int nChunk = (len + ARG.getIntArg("chunk") - 1) / ARG.getIntArg("chunk");
            for (int i = 0; i < nChunk; i++) {
                size_t s = i * len / nChunk, t = (i + 1) * len / nChunk;
                formatGene(ind2species, pos + s, t - s, len);
            }
        }
        return true;
    }

    Workflow(int argc, char** argv){
        //string mappingFile;
        meta.initialize(argc, argv);
        if (ARG.getStringArg("root") != "") addName(ARG.getStringArg("root"));
        init();
    }

    void init(){
        tripInit.nThreads = meta.nThreads;
		
		if (ARG.getIntArg("ambiguity") == 0){
			if (ARG.getStringArg("format") == "fasta") {
				readFasta(ARG.getStringArg("input"));
			}
			else if (ARG.getStringArg("format") == "list") {
				ifstream fin(ARG.getStringArg("input"));
				string line;
				while (getline(fin, line)) {
					readFasta(line);
				}
			}
			else if (ARG.getStringArg("format") == "phylip") {
				ifstream fin(ARG.getStringArg("input")), fin2(ARG.getStringArg("input"));
				string line;
				while (readPhylip(fin, fin2));
			}
			else {
				cerr << "Failed to parse format named '" << ARG.getStringArg("format") << "'\n";
				exit(1);
			}
		}
		else {
			if (ARG.getStringArg("format") == "fasta") {
				readFastaAmbiguity(ARG.getStringArg("input"));
			}
			else if (ARG.getStringArg("format") == "list") {
				ifstream fin(ARG.getStringArg("input"));
				string line;
				while (getline(fin, line)) {
					readFastaAmbiguity(line);
				}
			}
			else if (ARG.getStringArg("format") == "phylip") {
				ifstream fin(ARG.getStringArg("input")), fin2(ARG.getStringArg("input"));
				string line;
				while (readPhylip(fin, fin2));
			}
			else {
				cerr << "Failed to parse format named '" << ARG.getStringArg("format") << "'\n";
				exit(1);
			}
		}
        tripInit.nSpecies = names.size();
    }

    BinaryTree getGenetree(){
        LOG.enabled = false;
        pair<score_t, string> result = meta.run();
        LOG.enabled = true;
        return BinaryTree(result.second, name2id);
    }
};

int main(int argc, char** argv){
    ARG.setProgramName("caster-site", "Coalescence-aware Alignment-based Species Tree EstimatoR (Site)");
    ARG.addStringArg('f', "format", "fasta", "Input file type, fasta: one fasta file for the whole alignment, list: a txt file containing a list of FASTA files, phylip: a phylip file for the whole alignment", true);
    ARG.addStringArg('m', "mutation", "", "Substitution rate file from Iqtree if assumming heterogeneous rates", true);
    ARG.addIntArg('d', "diskcover", 1, "The number of replicates in the disk covering method", true);
	ARG.addIntArg(0, "ambiguity", 0, "0 (default): ambiguity codes are treated as N, 1: ambiguity codes are treated as diploid unphased sites");
    ARG.addIntArg(0, "chunk", 10000, "The chunk size of each local region for parameter estimation");
    
    Workflow WF(argc, argv);
    LOG << "#Base: " << WF.meta.tripInit.seq.len() << endl;
    auto res = WF.meta.run();
    LOG << "Score: " << (double) res.first << endl;
	return 0;
}
