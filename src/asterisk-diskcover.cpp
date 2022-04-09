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

#define ERROR_TOLERANCE 0.001

#include "argparser.hpp"
#include "sequtils.hpp"
#include "biallelic-diskcover.hpp"
#include "algorithms.hpp"

using namespace std;

const int NUM_BINS = 10, MIN_JURY_SIZE = 50, MAX_CANDIDATE_SIZE = 200, VIP_CANDIDATE_SIZE = 100, MIN_BIN_SIZE = 50;

MetaAlgorithm meta;
TripartitionInitializer &tripInit = meta.tripInit;

vector<string> &names = meta.names;
unordered_map<string, int> &name2id = meta.name2id;

double shortQuartetKeepRate = 0.95;
double rateModifier = 1;
int replicate = 3;

mutex mtx;

void addName(const string &name){
    if (name2id.count(name) == 0){
        name2id[name] = names.size();
        names.push_back(name);
    }
}

string formatName(const string &name){
	string res;
	for (char c: name){
		if (c != '>' && c != ' ' && c != '\t') res += c;
	}
	return res;
}

vector<pair<int, double> > infoSort(const vector<array<vector<bool>, 4> > &sites,
        const vector<int> &vips, const unordered_set<int> &singletonSet){
	mt19937 generator;
	vector<int> candidate, ordering, jury;
	unordered_set<int> vipset;
	for (int i: vips){
		ordering.push_back(i);
	}
	shuffle(ordering.begin(), ordering.end(), generator);
	for (int i = 0, imax = min(VIP_CANDIDATE_SIZE, (int) vips.size()); i < imax; i++){
		candidate.push_back(ordering[i]);
		vipset.insert(i);
	}
	ordering.clear();
	for (int i = 0; i < sites.size(); i++){
		if (vipset.count(i) == 0 && singletonSet.count(i) == 0) ordering.push_back(i);
	}
	shuffle(ordering.begin(), ordering.end(), generator);
	for (int i = 0, imax = min((int) ordering.size(), MAX_CANDIDATE_SIZE - (int) candidate.size()); i < imax; i++){
		candidate.push_back(ordering[i]);
	}
	vector<pair<int, double> > candidatePair, allPair;
	for (int i: candidate){
		double MI = 0;
		for (int j: candidate){
			if (i == j) continue;
			double freq[4][4] = {}, ifreq[4] = {}, jfreq[4] = {}, tfreq = 0;
			for (int x = 0; x < 4; x++){
				for (int y = 0; y < 4; y++){
					int cnt = 0;
					for (int p = 0; p < sites[i][x].size(); p++){
						cnt += (sites[i][x][p] & sites[j][y][p]);
					}
					freq[x][y] = cnt;
					ifreq[x] += freq[x][y];
					jfreq[y] += freq[x][y]; 
					tfreq += freq[x][y];
				}
			}
			for (int x = 0; x < 4; x++){
				for (int y = 0; y < 4; y++){
					if (freq[x][y] == 0) continue;
					MI += freq[x][y] * (log(freq[x][y] * tfreq / (ifreq[x] * jfreq[y])));
				}
			}
		}
		candidatePair.push_back({i, MI});
	}
	sort(candidatePair.begin(), candidatePair.end(), [](pair<int, double> a, pair<int, double> b){return a.second > b.second;});
	int jurySize = max(MIN_JURY_SIZE, MAX_CANDIDATE_SIZE * MAX_CANDIDATE_SIZE / (int) sites.size());
	for (int i = 0; i < min((int) candidate.size(), jurySize); i++){
		jury.push_back(candidatePair[i].first);
	}
	for (int i = 0; i < sites.size(); i++){
		double MI = 0;
		int juryCnt = 0;
		for (int j: jury){
			if (i == j) continue;
			juryCnt++;
			double freq[4][4] = {}, ifreq[4] = {}, jfreq[4] = {}, tfreq = 0;
			for (int x = 0; x < 4; x++){
				for (int y = 0; y < 4; y++){
					int cnt = 0;
					for (int p = 0; p < sites[i][x].size(); p++){
						cnt += (sites[i][x][p] & sites[j][y][p]);
					}
					freq[x][y] = cnt;
					ifreq[x] += freq[x][y];
					jfreq[y] += freq[x][y];
					tfreq += freq[x][y];
				}
			}
			for (int x = 0; x < 4; x++){
				for (int y = 0; y < 4; y++){
					if (freq[x][y] == 0) continue;
					MI += freq[x][y] / tfreq * log(freq[x][y] * tfreq / (ifreq[x] * jfreq[y]));
				}
			}
		}
		allPair.push_back({i, MI / juryCnt});
	}
	sort(allPair.begin(), allPair.end(), [](pair<int, double> a, pair<int, double> b){return a.second > b.second;});
	return allPair;
}

void addBin(const vector<int> &geneTaxa, const vector<vector<int> > &bins,
        const MSA &msa, const DistanceMatrix &guideDM, const Tree &guideTree){
	mt19937 randomEngine;
    TripartitionInitializer temp;
    int n = names.size();
    temp.taxa2block.resize(n);
    temp.taxa2row.resize(n);
	for (const vector<int> &bin: bins){
        const MSA msaBin = msa.sample(bin);
        //cerr << (string) msaBin;
        //cerr << msaBin.nSites() << endl;
        if (msaBin.nSites() < 5) continue;
        AlignmentHot a[7], b[7];
        a[0] = msaBin.A(); b[0] = msaBin.C() + msaBin.G() + msaBin.T();
        a[1] = msaBin.C(); b[1] = msaBin.A() + msaBin.G() + msaBin.T();
        a[2] = msaBin.G(); b[2] = msaBin.A() + msaBin.C() + msaBin.T();
        a[3] = msaBin.T(); b[3] = msaBin.A() + msaBin.C() + msaBin.G();
        a[4] = msaBin.A() + msaBin.C(); b[4] = msaBin.G() + msaBin.T();
        a[5] = msaBin.A() + msaBin.G(); b[5] = msaBin.C() + msaBin.T();
        a[6] = msaBin.A() + msaBin.T(); b[6] = msaBin.C() + msaBin.G();
        for (int r = 0; r < 7; r++){
            DistanceMatrix dm(a[r], b[r], false);
            double sumGuideDM = 0, sumDM = 0;
            for (int i = 0; i < n; i++){
                for (int j = 0; j < n; j++){
                    if (dm[i][j] < 2 && guideDM[i][j] < 2){
                        sumDM += dm[i][j];
                        sumGuideDM += guideDM[i][j];
                    }
                }
            }
            if (sumGuideDM <= 1e-6) continue;
            const vector<vector<int> > dc = guideTree.diskCovering(randomEngine, rateModifier * sumDM / sumGuideDM, replicate);
            
            for (const vector<int> &group: dc){
                if (group.size() < 4) continue;
                AlignmentHot block1, block2;
                block1.p = a[r].p; block2.p = b[r].p;
                for (int i = 0; i < group.size(); i++){
                    int j = group[i];
                    block1.s.push_back(a[r][j]);
                    block2.s.push_back(b[r][j]);
                    temp.taxa2block[geneTaxa[j]].push_back(temp.block1.size());
                    temp.taxa2row[geneTaxa[j]].push_back(i);
                }
                temp.block1.push_back(block1);
                temp.block2.push_back(block2);
                temp.w.push_back(1.0 / (msaBin.nSites() - 1));
            }
        }
    }
    {
	    const lock_guard<mutex> lock(mtx);
	    int offset = tripInit.w.size();
        for (int i = 0; i < temp.w.size(); i++){
            tripInit.block1.push_back(temp.block1[i]);
            tripInit.block2.push_back(temp.block2[i]);
            tripInit.w.push_back(temp.w[i]);
        }
        for (int i = 0; i < temp.taxa2block.size(); i++){
            for (int j: temp.taxa2block[i]) tripInit.taxa2block[i].push_back(j + offset);
            for (int j: temp.taxa2row[i]) tripInit.taxa2row[i].push_back(j);
        }
    }
}


void formatGene(const vector<int> geneTaxa, const vector<string> geneSeqs){
	MSA msa(geneSeqs);
    DistanceMatrix guideDM(msa.A()+msa.G(), msa.C()+msa.T());
    Tree guideTree(guideDM);
    
	vector<int> triallelic, alleleCnt = msa.alleleCnt(), minorAlleleSum = msa.minorAlleleSum();
    unordered_set<int> singletonSet, biallelicSet, triallelicSet, quadriallelicSet;
    for (int i = 0; i < alleleCnt.size(); i++){
        if (alleleCnt[i] == 3) triallelic.push_back(i);
        if (minorAlleleSum[i] <= 1) singletonSet.insert(i);
        if (alleleCnt[i] == 2 && minorAlleleSum[i] > 1) biallelicSet.insert(i);
        if (alleleCnt[i] == 3) triallelicSet.insert(i);
        if (alleleCnt[i] == 4) quadriallelicSet.insert(i);
    }
    vector<array<vector<bool>, 4> > sites(msa.nSites());
	for (int i = 0; i < msa.nSites(); i++){
        for (int x = 0; x < 4; x++){
            for (int j = 0; j < msa.nTaxa(); j++){
                sites[i][x].push_back(msa[x][j][i]);
            }
        }
    }
    vector<pair<int, double> > ordered = infoSort(sites, triallelic, singletonSet);
    vector<int> orderedDivision[2];
    for (const pair<int, double> &e: ordered){
        if (biallelicSet.count(e.first) || triallelicSet.count(e.first)) orderedDivision[0].push_back(e.first);
        if (quadriallelicSet.count(e.first)) orderedDivision[1].push_back(e.first);
    }
    vector<vector<int> > bins;
    for (vector<int> &division: orderedDivision){
        int nBin = max(1, min(NUM_BINS, (int) division.size() / MIN_BIN_SIZE));
        for (int i = 0; i < nBin; i++){
            vector<int> bin;
            for (int j = i * division.size() / nBin, jmax = (i + 1) * division.size() / nBin; j < jmax; j++){
                bin.push_back(division[j]);
            }
            bins.push_back(bin);
        }
    }
    addBin(geneTaxa, bins, msa, guideDM, guideTree);
}
/*
void addBin(const vector<int> &geneTaxa, const vector<vector<int> > &bins, const MSA &msa){
	mt19937 randomEngine;
    TripartitionInitializer temp;
    int n = names.size();
    temp.taxa2block.resize(n);
    temp.taxa2row.resize(n);
	for (const vector<int> &bin: bins){
        const MSA msaBin = msa.sample(bin);
        AlignmentHot a[7], b[7];
        a[0] = msaBin.A(); b[0] = msaBin.C() + msaBin.G() + msaBin.T();
        a[1] = msaBin.C(); b[1] = msaBin.A() + msaBin.G() + msaBin.T();
        a[2] = msaBin.G(); b[2] = msaBin.A() + msaBin.C() + msaBin.T();
        a[3] = msaBin.T(); b[3] = msaBin.A() + msaBin.C() + msaBin.G();
        a[4] = msaBin.A() + msaBin.C(); b[4] = msaBin.G() + msaBin.T();
        a[5] = msaBin.A() + msaBin.G(); b[5] = msaBin.C() + msaBin.T();
        a[6] = msaBin.A() + msaBin.T(); b[6] = msaBin.C() + msaBin.G();
        for (int r = 0; r < 7; r++){
            double rate = 0; // log((double)n)/n
            const vector<vector<int> > dc = diskCovering(a[r], b[r], randomEngine, rate, 2, (int) (2*exp(2.5*rate)));
            for (const vector<int> &group: dc){
                if (group.size() < 4) continue;
                AlignmentHot block1, block2;
                block1.p = a[r].p; block2.p = b[r].p;
                for (int i = 0; i < group.size(); i++){
                    int j = group[i];
                    block1.s.push_back(a[r][j]);
                    block2.s.push_back(b[r][j]);
                    temp.taxa2block[geneTaxa[j]].push_back(temp.block1.size());
                    temp.taxa2row[geneTaxa[j]].push_back(i);
                }
                temp.block1.push_back(block1);
                temp.block2.push_back(block2);
                temp.w.push_back(1);
            }
        }
    }
    {
	    const lock_guard<mutex> lock(mtx);
	    int offset = tripInit.w.size();
        for (int i = 0; i < temp.w.size(); i++){
            tripInit.block1.push_back(temp.block1[i]);
            tripInit.block2.push_back(temp.block2[i]);
            tripInit.w.push_back(temp.w[i]);
        }
        for (int i = 0; i < temp.taxa2block.size(); i++){
            for (int j: temp.taxa2block[i]) tripInit.taxa2block[i].push_back(j + offset);
            for (int j: temp.taxa2row[i]) tripInit.taxa2row[i].push_back(j);
        }
    }
}

void formatGene(const vector<int> geneTaxa, const vector<string> geneSeqs){
	MSA msa(geneSeqs);
    DistanceMatrix DM(msa.A()+msa.G(), msa.C()+msa.T(), false);
    vector<int> nameInt;
    for (string nm: names) nameInt.push_back(stoi(nm));
    Tree T(DM, nameInt);
    cerr << (string) T;
    exit(0); 

	vector<int> triallelic, alleleCnt = msa.alleleCnt(), minorAlleleSum = msa.minorAlleleSum();
    unordered_set<int> singletonSet, biallelicSet, triallelicSet, quadriallelicSet;
    for (int i = 0; i < alleleCnt.size(); i++){
        if (alleleCnt[i] == 3) triallelic.push_back(i);
        if (minorAlleleSum[i] <= 1) singletonSet.insert(i);
        if (alleleCnt[i] == 2 && minorAlleleSum[i] > 1) biallelicSet.insert(i);
        if (alleleCnt[i] == 3) triallelicSet.insert(i);
        if (alleleCnt[i] == 4) quadriallelicSet.insert(i);
    }
    vector<array<vector<bool>, 4> > sites(msa.nSites());
	for (int i = 0; i < msa.nSites(); i++){
        for (int x = 0; x < 4; x++){
            for (int j = 0; j < msa.nTaxa(); j++){
                sites[i][x].push_back(msa[x][j][i]);
            }
        }
    }
    vector<pair<int, double> > ordered = infoSort(sites, triallelic, singletonSet);
    int nBin = max(1, min(NUM_BINS, (int) ordered.size() / MIN_BIN_SIZE));
    vector<int> orderedDivision[2];
    for (const pair<int, double> &e: ordered){
        if (biallelicSet.count(e.first) || triallelicSet.count(e.first)) orderedDivision[0].push_back(e.first);
        if (quadriallelicSet.count(e.first)) orderedDivision[1].push_back(e.first);
    }
    vector<vector<int> > bins;
    for (vector<int> &division: orderedDivision){
        for (int i = 0; i < nBin; i++){
            vector<int> bin;
            for (int j = i * division.size() / nBin, jmax = (i + 1) * division.size() / nBin; j < jmax; j++){
                bin.push_back(division[j]);
            }
            bins.push_back(bin);
        }
    }
    addBin(geneTaxa, bins, msa);
}
*/
void preprocessPhilip(istream &fin){
	int M, L, n = 0;
	vector<thread> thrds;
    vector<double> tempThreshold(meta.nThreads), threshold;

	while (fin >> M){
		fin >> L;
        vector<string> geneSeqs;
		for (int i = 0; i < M; i++){
			string s, line;
			fin >> s >> line;
			geneSeqs.push_back(line);
		}
        thrds.emplace_back([](double &result, const vector<string> seqs, double keepRate){
            MSA msa(seqs);
            DistanceMatrix guideDM(msa.A()+msa.G(), msa.C()+msa.T());
            Tree guideTree(guideDM);
            vector<double> lengths = guideTree.shortQuartetLengths();
            sort(lengths.begin(), lengths.end());
            result = lengths[int((lengths.size()-1) * keepRate)];
        }, ref(tempThreshold[n]), move(geneSeqs), shortQuartetKeepRate);
		n++;
		if (n == meta.nThreads){
			for (thread &t: thrds) t.join();
            for (int i = 0; i < n; i++) threshold.push_back(tempThreshold[i]);
			thrds.clear();
            n = 0;
		}
	}
	for (thread &t: thrds) t.join();
    for (int i = 0; i < n; i++) threshold.push_back(tempThreshold[i]);
    sort(threshold.begin(), threshold.end());
    rateModifier = log(2) / threshold[threshold.size() / 2];
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
		thrds.emplace_back(formatGene, move(geneTaxa), move(geneSeqs));
		if (thrds.size() == meta.nThreads){
			for (thread &t: thrds) t.join();
			thrds.clear();
		}
	}
	for (thread &t: thrds) t.join();
}

void preprocessFasta(string file, int &cnt, double &qsum){
	ifstream fin(file);
	vector<string> geneSeqs;
	string line;
	while (getline(fin, line)){
		if (line[0] == '>') geneSeqs.emplace_back();
		else geneSeqs.back() += formatName(line);
	}
	MSA msa(geneSeqs);
    DistanceMatrix guideDM(msa.A()+msa.G(), msa.C()+msa.T());
    Tree guideTree(guideDM);
    for (double v: guideTree.shortQuartetLengths()){
        cnt++;
        qsum += v;
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
	ARG.setProgramName("asterisk-diskcover", "Accurate Species Tree EstimatoR from Independent Site Kernals");
	ARG.addStringArg('f', "format", "fasta", "Input file type, fasta: a txt file containing a list of FASTA gene files; phylip: a phylip file with all genes");
	ARG.addFlag('F', "formatphylip", "Indicating input file as a PHYLIP file containing all genes (`-f phylip`)", [&](){
		ARG.getStringArg("format") = "phylip";
	}, true);
	string mappingFile;
	meta.initialize(argc, argv, " -y", HELP_TEXT);
	tripInit.numThreads = meta.nThreads;

	bool phylip = ARG.getStringArg("format").compare("fasta");
	
	if (phylip) {
        int M, L;
        {
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
        tripInit.taxa2block.resize(names.size());
        tripInit.taxa2row.resize(names.size());
        {
		    ifstream fin(ARG.getStringArg("input"));
		    preprocessPhilip(fin);
        }
        {
		    ifstream fin(ARG.getStringArg("input"));
		    readPhilip(fin);
        }
	}
	else {
		ifstream listIn(ARG.getStringArg("input"));
        vector<string> fileNames;
		vector<thread> thrds;
		for (string file; getline(listIn, file);){
            fileNames.push_back(file);
        }
        for (string file: fileNames){
            ifstream fin(file);
            string line;
            while (getline(fin, line)){
                if (line[0] == '>') addName(formatName(line));
            }
		}
        tripInit.taxa2block.resize(names.size());
        tripInit.taxa2row.resize(names.size());
        vector<int> cnt(meta.nThreads);
        vector<double> qSum(meta.nThreads);
        int n = 0;
		for (string file: fileNames){
			thrds.emplace_back(preprocessFasta, file, ref(cnt[n]), ref(qSum[n]));
			n++;
            if (n == meta.nThreads){
				for (thread &t: thrds) t.join();
				thrds.clear();
			}
		}
        int totalCnt = 0;
        double totalqSum = 0;
        for (int v: cnt) totalCnt += v;
        for (double v: qSum) totalqSum += v;
        rateModifier = log(2) * totalCnt / totalqSum;
        for (string file: fileNames){
			thrds.emplace_back(readFasta, file);
			if (thrds.size() == meta.nThreads){
				for (thread &t: thrds) t.join();
				thrds.clear();
			}
		}
		for (thread &t: thrds) t.join();
	}
	
    cerr << "rateModifier: " << rateModifier << endl;
	cerr << "#Blocks: " << tripInit.block1.size() << endl;
	
	score_t score = meta.run().first;
	cerr << "Score: " << (double) score << endl;
	return 0;
}
