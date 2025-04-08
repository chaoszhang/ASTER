#ifndef SEQ_UTILS
#define SEQ_UTILS

#include<queue>
#include<tuple>
#include<algorithm>
#include<random>
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<unordered_map>
#include<unordered_set>

using namespace std;

namespace SeqUtils{
    string QUALITY2ASCII = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";

    string fastaFormatName(const string &name){
        string res;
        for (char c: name){
            if (c != '>' && c != ' ' && c != '\t') res += c;
        }
        return res;
    }

    string fastaFormatRead(const string &name){
        string res;
        for (char c: name){
            if (('a' <= c && 'z' >= c) || ('A' <= c && 'Z' >= c) || c == '-') res += c;
        }
        return res;
    }
}

struct SeqParser{
    bool isFasta = false;
    string seq, quality, L1;
    ifstream fin;

    const long long STEP_SIZE = 1 << 27;
    long long cnt = 0, lenCnt = 0, threshold = STEP_SIZE;
    
    SeqParser(string fileName, string fileFormat = "auto"){
        ifstream ftemp(fileName);
        string temp;
        if (!getline(ftemp, temp)){
            cerr << "Input file path '" << fileName << "' seems empty!\n";
            exit(0);
        }
        if (temp.length() == 0){
            cerr << "Input file path '" << fileName << "' seems empty!\n";
            exit(0);
        }
        if (fileFormat == "fastq" || (fileFormat == "auto" && seemsFastq(temp))) initFastq(fileName);
        else if (fileFormat == "fasta" || (fileFormat == "auto" && seemsFasta(temp))) initFasta(fileName);
        else {
            cerr << "Input file '" << fileName << "' bad format! Only FASTA and FASTQ are supported! (Did you forget unzipping it?)\n";
            exit(0);
        }
    }

    bool nextSeq(){
        if (isFasta) return nextFastaSeq();
        else return nextFastqSeq();
    }

    string getSeq(size_t threshold = 0){
        if (isFasta) return seq;
        if (threshold >= SeqUtils::QUALITY2ASCII.length()) threshold = SeqUtils::QUALITY2ASCII.length() - 1;
        string res;
        for (size_t i = 0; i < seq.size(); i++){
            if (quality[i] >= SeqUtils::QUALITY2ASCII[threshold]) res += seq[i];
            else res += 'N';
        }
        return res;
    }

    static string STRING_REVERSE_COMPLEMENT(const string &seq){
        string res;
        for (int i = seq.size() - 1; i >= 0; i--){
            switch (seq[i]) {
                case 'A': res += 'T'; break;
                case 'C': res += 'G'; break;
                case 'G': res += 'C'; break;
                case 'T': res += 'A'; break;
                default: res += 'N';
            }
        }
        return res;
    }

    static string formatSeq(const string raw){
        string res;
        for (char c: raw){
            switch (c) {
                case 'A': case 'a': res += 'A'; break;
                case 'C': case 'c': res += 'C'; break;
                case 'G': case 'g': res += 'G'; break;
                case 'T': case 't': case 'U': case 'u': res += 'T'; break;
                case ' ': case '\t': case '\r': case '\n': break; 
                default: res += 'N';
            }
        }
        return res;
    }

private:
    static bool seemsFastq(const string &line){
        return line[0] == '@';
    }

    static bool seemsFasta(const string &line){
        return line[0] == '>';
    }

    void initFastq(const string &fileName){
        isFasta = false;
        fin.open(fileName);
    }

    void initFasta(const string &fileName){
        isFasta = true;
        fin.open(fileName);
        if (!getline(fin, L1)){
            cerr << "Input file path '" << fileName << "' seems empty!\n";
            exit(0);
        }
        if (L1.length() == 0) getline(fin, L1);
    }

    bool nextFastaSeq(){
        if (L1 == "") return false;

        string line;
        seq = "";
        while (getline(fin, line)){
            if (line.length() > 0 && line[0] == '>') {
                cnt++;
                lenCnt += seq.length();
                if (lenCnt >= threshold){
                    while (lenCnt >= threshold) threshold += STEP_SIZE;
                    LOG << cnt << " sequences read (" << lenCnt << " BPs); the last sequence is '" + L1 + "'" << endl;
                }

                L1 = line;
                return true;
            }
            seq += formatSeq(line);
        }
        
        cnt++;
        lenCnt += seq.length();
        if (lenCnt >= threshold){
            while (lenCnt >= threshold) threshold += STEP_SIZE;
            LOG << cnt << " sequences read (" << lenCnt << " BPs); the last sequence is '" + L1 + "'" << endl;
        }

        L1 = "";
        return true;
    }

    bool nextFastqSeq(){
        string L2, L3;
        if (!getline(fin, L1)) return false;
        if (L1.length() == 0) getline(fin, L1);
        getline(fin, L2);
        getline(fin, L3);
        getline(fin, quality);
        seq = formatSeq(L2);

        cnt++;
        lenCnt += seq.length();
        if (lenCnt >= threshold){
            while (lenCnt >= threshold) threshold += STEP_SIZE;
            LOG << cnt << " sequences read (" << lenCnt << " BPs); the last sequence is '" + L1 + "'" << endl;
        }
        return true;
    }
};

struct AlignmentParser{
    bool isFasta = false, isFastaList = false, isPhylip = false;
    bool formatAmbiguity = false, formatAA2NA = false, formatAA = false, formatNA = false;

    queue<string> fastaQueue;
    string bufferName, bufferSeq, nextLine;
    size_t length;
    int line, phylipNspecies;
    bool firstFastaSeq;
    ifstream fin;

    bool ambiguitySecond = true;
    
    AlignmentParser(string fileName, string fileFormat = "auto", string seqFormat = "NA"){
        ifstream ftemp(fileName);
        string temp;
        if (!getline(ftemp, temp)){
            cerr << "Input file path '" << fileName << "' seems empty!\n";
            exit(0);
        }
        if (temp.length() == 0){
            cerr << "Input file path '" << fileName << "' seems empty!\n";
            exit(0);
        }
        if (fileFormat == "phylip" || (fileFormat == "auto" && seemsPhylip(temp))) initPhylip(fileName);
        else if (fileFormat == "fasta" || (fileFormat == "auto" && seemsFasta(temp))) initFasta(fileName);
        else initFastaList(fileName);
        if (seqFormat == "ambiguity") formatAmbiguity = true;
        else if (seqFormat == "AA2NA") formatAA2NA = true;
        else if (seqFormat == "AA") formatAA = true;
        else formatNA = true;
    }

    bool nextAlignment(){
        ambiguitySecond = true;
        if (isFasta) return nextAlignmentFasta();
        else if (isPhylip) return nextAlignmentPhylip();
        else return false;
    }

    int getLength(){
        return length;
    }

    bool nextSeq(){
        if (formatAmbiguity){
            ambiguitySecond = (!ambiguitySecond);
            if (ambiguitySecond) return true;
        }

        if (isFasta) return nextSeqFasta();
        else if (isPhylip) return nextSeqPhylip();
        else return false;
    }

    string getName(){
        return bufferName;
    }

    string getSeq(){
        if (bufferSeq.length() != length){
            cerr << "Error: Sequence length of " << bufferName << " is not consistent!\n";
            exit(0);
        }
        if (formatNA) return toFormatNA(bufferSeq);
        if (formatAmbiguity) return toAmbiguity(bufferSeq);
        if (formatAA2NA) return toFormatAA2NA(bufferSeq);
        return bufferSeq;
    }

private:
    static string getFastaName(const string &fasta){
        size_t i = 0;
        string res;
        // while (i < fasta.size() && fasta[i] == ' ') i++;
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

    static string removeSpace(const string &seq){
        string res;
        for (char c: seq){
            if (c != ' ' && c != '\t' && c != '\r') res += c;
        }
        return res;
    }

    static bool seemsPhylip(const string &line){
        size_t i = 0;
        while (i < line.length() && (line[i] == ' ' || line[i] == '\t')) i++;
        if (i == line.length() || line[i] < '0' || line[i] > '9') return false;
        while (i < line.length() && line[i] >= '0' && line[i] <= '9') i++;
        if (i == line.length() || (line[i] != ' ' && line[i] != '\t')) return false;
        while (i < line.length() && (line[i] == ' ' || line[i] == '\t')) i++;
        if (i == line.length() || line[i] < '0' || line[i] > '9') return false;
        while (i < line.length() && line[i] >= '0' && line[i] <= '9') i++;
        while (i < line.length() && (line[i] == ' ' || line[i] == '\t')) i++;
        return i == line.length();
    }

    static bool seemsFasta(const string &line){
        return line[0] == '>';
    }

    void initPhylip(const string &fileName){
        isPhylip = true;
        fin.open(fileName);
    }

    void initFasta(const string &fileName){
        isFasta = true;
        fastaQueue.push(fileName);
    }

    void initFastaList(const string &fileName){
        isFasta = true;
        ifstream ftemp(fileName);
        string temp;
        while (getline(ftemp, temp)){
            fastaQueue.push(temp);
        }
    }

    bool parseSeqFasta(){
        bufferName = "";
        bufferSeq = "";
        string line;
        if (nextLine.length()) {
            line = nextLine;
            nextLine = "";
        }
        else if (!getline(fin, line)) return false;

        bufferName = getFastaName(line);
        while (getline(fin, line)){
            if (line[0] == '>') {
                nextLine = line;
                break;
            }
            bufferSeq += removeSpace(line);
        }
        return true;
    }

    bool nextAlignmentFasta(){
        if (fastaQueue.empty()) return false;
        string fileName = fastaQueue.front();
        fastaQueue.pop();
        cerr << "Processing " << fileName << " ... \n";
        fin.close();
        fin.open(fileName);
        if (!parseSeqFasta()) {
            cerr << "Error: FASTA file empty!\n";
            exit(0);
        }
        length = bufferSeq.size();
        firstFastaSeq = true;
        return true;
    }

    bool nextSeqFasta(){
        if (firstFastaSeq) {
            firstFastaSeq = false;
            return true;
        }
        return parseSeqFasta();
    }

    bool nextAlignmentPhylip(){
        if (!(fin >> phylipNspecies)) return false;
        fin >> length;
        cerr << "Reading alignment of " << phylipNspecies << " species of length " << length << " ...\n";
        return true;
    }

    bool nextSeqPhylip(){
        if (phylipNspecies == 0) return false;
        phylipNspecies--;
        fin >> bufferName >> bufferSeq;
        bufferSeq = removeSpace(bufferSeq);
        return true;
    }

    string toFormatNA(const string &seq){
        string res;
        for (char c: seq){
            switch (c) {
                case 'A': case 'a': res += 'A'; break;
                case 'C': case 'c': res += 'C'; break;
                case 'G': case 'g': res += 'G'; break;
                case 'T': case 't': case 'U': case 'u': res += 'T'; break;
                default: res += '-';
            }
        }
        return res;
    }

    string toAmbiguity(const string &seq){
        string res;
        if (!ambiguitySecond){
            for (char c: seq){
                switch (c) {
                    case 'A': case 'a': case 'M': case 'm': case 'R': case 'r': case 'W': case 'w': res += 'A'; break;
                    case 'C': case 'c': case 'S': case 's': case 'Y': case 'y': res += 'C'; break;
                    case 'G': case 'g': case 'K': case 'k': res += 'G'; break;
                    case 'T': case 't': case 'U': case 'u': res += 'T'; break;
                    default: res += '-';
                }
            }
        }
        else {
            for (char c: seq){
                switch (c) {
                    case 'A': case 'a': res += 'A'; break;
                    case 'C': case 'c': case 'M': case 'm': res += 'C'; break;
                    case 'G': case 'g': case 'R': case 'r': case 'S': case 's': res += 'G'; break;
                    case 'T': case 't': case 'U': case 'u': case 'W': case 'w': case 'Y': case 'y': case 'K': case 'k': res += 'T'; break;
                    default: res += '-';
                }
            }
        }
        return res;
    }

    string toFormatAA2NA(const string &seq){
        string res;
        for (char c: seq){
            switch (c) {
                case 'C': case 'c': case 'M': case 'm': case 'I': case 'i': case 'L': case 'l': case 'V': case 'v': res += 'A'; break;
			    case 'D': case 'd': case 'E': case 'e': case 'Q': case 'q': case 'N': case 'n': case 'H': case 'h': case 'R': case 'r': case 'K': case 'k': res += 'T'; break;
			    case 'S': case 's': case 'T': case 't': case 'A': case 'a': case 'G': case 'g': case 'P': case 'p': res += 'C'; break;
			    case 'W': case 'w': case 'Y': case 'y': case 'F': case 'f': res += 'G'; break;
                default: res += '-';
            }
        }
        return res;
    }
};

struct SeqHot{
    vector<bool> b;

    SeqHot(){}
    SeqHot(int n): b(n){}
    SeqHot(const vector<bool> b): b(b){}

    bool operator[](int i) const{
        return b[i];
    }

    int nSites() const{
        return b.size(); 
    }

    SeqHot operator|(const SeqHot x) const{
        int n = b.size();
        SeqHot y(n);
        for (int i = 0; i < n; i++){
            y.b[i] = (b[i] | x.b[i]);
        }
        return y;
    }

    SeqHot operator&(const SeqHot x) const{
        int n = b.size();
        SeqHot y(n);
        for (int i = 0; i < n; i++){
            y.b[i] = (b[i] & x.b[i]);
        }
        return y;
    }

    SeqHot operator^(const SeqHot x) const{
        int n = b.size();
        SeqHot y(n);
        for (int i = 0; i < n; i++){
            y.b[i] = (b[i] ^ x.b[i]);
        }
        return y;
    }

    operator string() const{
        string s;
        for (bool e: b) s += to_string(e);
        s += "\n";
        return s;
    }

    SeqHot sample(const vector<int> &cols) const{
        SeqHot result;
        for (int col: cols){
            result.b.push_back(b[col]);
        }
        return result;
    }

    int sum() const{
        int result = 0;
        for (bool v: b) result += v;
        return result;
    }
};

struct AlignmentHot
{
    vector<SeqHot> s;
    double p;

    AlignmentHot(){}
    AlignmentHot(int n): s(n){}
    AlignmentHot(const vector<SeqHot> s, double p): s(s), p(p){}

    SeqHot& operator[](int i){
        return s[i];
    }

    const SeqHot& operator[](int i) const{
        return s[i];
    }

    int nTaxa() const{
        return s.size(); 
    }

    int nSites() const{
        return s[0].nSites();
    }

    AlignmentHot operator+(const AlignmentHot &x) const{
        AlignmentHot y;
        for (size_t i = 0; i < s.size(); i++){
            y.s.push_back(s[i] | x[i]);
        }
        y.p = p + x.p;
        return y;
    }

    AlignmentHot operator-(const AlignmentHot &x) const{
        AlignmentHot y;
        for (size_t i = 0; i < s.size(); i++){
            y.s.push_back(s[i] ^ x[i]);
        }
        y.p = p - x.p;
        return y;
    }

    operator string() const{
        string str = "p = ";
        str += to_string(p) + "\n";
        for (const SeqHot &e: s) str += (string) e;
        return str;
    }

    SeqHot seqOr() const{
        SeqHot result = s[0];
        for (size_t i = 1; i < s.size(); i++)
            result = result | s[i];
        return result;
    }

    vector<int> seqSum() const{
        vector<int> result(nSites());
        for (int i = 0; i < nTaxa(); i++)
            for (int j = 0; j < nSites(); j++)
                result[j] += s[i][j];
        return result;
    }

    AlignmentHot sample(const vector<int> &cols) const{
        AlignmentHot result;
        result.p = p;
        for (int i = 0; i < nTaxa(); i++){
            result.s.push_back(s[i].sample(cols));
        }
        return result;
    }
};

struct MSA
{
    AlignmentHot a[4];
    
    AlignmentHot& A(){ return a[0]; }
    AlignmentHot& C(){ return a[1]; }
    AlignmentHot& G(){ return a[2]; }
    AlignmentHot& T(){ return a[3]; }
    const AlignmentHot& A() const{ return a[0]; }
    const AlignmentHot& C() const{ return a[1]; }
    const AlignmentHot& G() const{ return a[2]; }
    const AlignmentHot& T() const{ return a[3]; }
    AlignmentHot& operator[](int i){ return a[i]; }
    const AlignmentHot& operator[](int i) const{ return a[i]; }

    MSA(){}

    MSA(const vector<string> seqs){
        int n = seqs.size(), L = seqs[0].size();
        int cntA = 0, cntC = 0, cntG = 0, cntT = 0; 
        for (int r = 0; r < 4; r++){
            for (int i = 0; i < n; i++){
                a[r].s.emplace_back(L);
            }
        }
        for (int i = 0; i < n; i++){
            for (int j = 0; j < L; j++){
                switch (seqs[i][j]){
                    case 'A': case 'a': cntA++; A().s[i].b[j] = true; break; 
                    case 'C': case 'c': cntC++; C().s[i].b[j] = true; break;
                    case 'G': case 'g': cntG++; G().s[i].b[j] = true; break;
                    case 'T': case 't': cntT++; T().s[i].b[j] = true; break;
                }
            }
        }
        double total = cntA + cntC + cntG + cntT; 
        A().p = cntA / total; C().p = cntC / total; G().p = cntG / total; T().p = cntT / total;
    }

    int nTaxa() const{
        return a[0].nTaxa(); 
    }

    int nSites() const{
        return a[0].nSites();
    }

    operator string() const{
        string str = "A: ";
        str += (string) A();
        str += "C: ";
        str += (string) C();
        str += "G: ";
        str += (string) G();
        str += "T: ";
        str += (string) T();
        return str;
    }
    
    vector<int> alleleCnt() const{
        SeqHot seqA = A().seqOr(), seqC = C().seqOr(), seqG = G().seqOr(), seqT = T().seqOr();
        vector<int> result;
        for (int i = 0; i < seqA.nSites(); i++){
            result.push_back(seqA[i] + seqC[i] + seqG[i] + seqT[i]);
        }
        return result;
    }

    vector<int> minorAlleleSum() const{
        vector<int> result(nSites()), maxS(nSites());
        for (int r = 0; r < 4; r++){
            vector<int> s = a[r].seqSum();
            for (int j = 0; j < nSites(); j++){
                result[j] += s[j];
                maxS[j] = max(maxS[j], s[j]);
            }
        }
        for (int j = 0; j < nSites(); j++){
            result[j] -= maxS[j];
        }
        return result;
    }

    MSA sample(const vector<int> &cols) const{
        MSA result;
        for (int i = 0; i < 4; i++){
            result[i] = a[i].sample(cols);
        }
        return result;
    }
};

struct DistanceMatrix{
    vector<vector<double> > d, w;
    
    vector<double>& operator[](int i){
        return d[i];
    }

    const vector<double>& operator[](int i) const{
        return d[i];
    }

    int nTaxa() const{
        return d.size(); 
    }

    DistanceMatrix(){}

    DistanceMatrix(int n): d(n, vector<double>(n)), w(n, vector<double>(n)){}

    DistanceMatrix(const AlignmentHot &a, const AlignmentHot &b, bool floyd = true): DistanceMatrix(a.nTaxa()){
        int n = a.nTaxa();
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                int x = ((a[i] & b[j]) | (a[j] & b[i])).sum();
                int y = ((a[i] | b[i]) & (a[j] | b[j])).sum();
                w[i][j] = y - x / (2 * a.p * b.p);
                if (w[i][j] > 0) d[i][j] = -log(w[i][j] / y);
                else {
                    w[i][j] = 0;
                    d[i][j] = 10;
                }
            }
        }
        if (!floyd) return;
        for (int k = 0; k < n; k++){
            for (int i = 0; i < n; i++){
                for (int j = 0; j < n; j++){
                    w[i][j] = min(w[i][j], w[i][k] + w[k][j]);
                }
            }
        }
    }

    vector<vector<int> > cluster(mt19937 &eng, double threshold = 3, bool isMax = true){
        vector<vector<int> > result;
        uniform_int_distribution<int> distribution(0,1);
        int n = nTaxa();
        vector<int> clusterNum;
        vector<tuple<double, int, int> > order;
        int nClusters = n, p = 0;
        for (int i = 0; i < n; i++){
            clusterNum.push_back(i);
        }
        for (int i = 0; i < n; i++){
            for (int j = i + 1; j < n; j++){
                order.push_back({(double) d[i][j], i, j});
            }
        }
        sort(order.begin(), order.end());
        while (nClusters > 1){
            int i, j;
            do { i = clusterNum[get<1>(order[p])]; j = clusterNum[get<2>(order[p])]; p++; } 
            while (i == -1 || j == -1 || i == j);
            vector<int> Si, Sj;
            double maxd = 0, mind = 10;
            for (int k = 0; k < n; k++){
                if (clusterNum[k] == i) Si.push_back(k);
                if (clusterNum[k] == j) Sj.push_back(k);
            }
            for (int ei: Si){
                for (int ej: Sj){
                    maxd = max(d[ei][ej], maxd);
                    mind = min(d[ei][ej], mind);
                }
            }
            if ((isMax && maxd > threshold) || (!isMax && mind > threshold)){
                if (distribution(eng)){
                    result.push_back(Si);
                    for (int k = 0; k < n; k++){
                        if (clusterNum[k] == i) clusterNum[k] = -1;
                    }
                }
                else {
                    result.push_back(Sj);
                    for (int k = 0; k < n; k++){
                        if (clusterNum[k] == j) clusterNum[k] = -1;
                    }
                }
                nClusters--;
            }
            else {
                for (int k = 0; k < n; k++){
                    if (clusterNum[k] == i) clusterNum[k] = j;
                }
                nClusters--;
            }
        }
        if (nClusters == 1){
            vector<int> S;
            for (int k = 0; k < n; k++){
                if (clusterNum[k] != -1) S.push_back(k);
            }
            result.push_back(S);
        }
        return result;
    }

    vector<DistanceMatrix> partition(const vector<vector<int> > &groups){
        vector<DistanceMatrix> result;
        for (const vector<int> &arr: groups){
            int n = arr.size();
            DistanceMatrix cur(n);
            for (int i = 0; i < n; i++){
                for (int j = 0; j < n; j++){
                    cur.d[i][j] = d[arr[i]][arr[j]];
                    cur.w[i][j] = w[arr[i]][arr[j]];
                }
            }
            result.push_back(cur);
        }
        return result;
    }

    operator string() const{
        string str;
        int n = nTaxa();
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                if (j != 0) str += "\t";
                str += to_string(d[i][j]);
            }
            str += "\n";
        }
        return str;
    }
};

#endif