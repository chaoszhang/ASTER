#ifndef SEQ_UTILS
#define SEQ_UTILS

#include<queue>
#include<tuple>
#include<algorithm>
#include<random>
#include<string>
#include<fstream>
#include<sstream>

namespace SeqUtils{
    string PARSE_LEAFNAME(const string& TEXT){
        string s;
        for (int i = 0; i < TEXT.size(); i++){
            if (TEXT[i] != '\"' && TEXT[i] != '\'') s += TEXT[i];
        }
        return s;
    }

    long long choose4(long long n){
        return n * (n - 1) * (n - 2) * (n - 3) / 24;
    }

    score_t from_string(const string s){
    	return stod(s);
    }

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

    vector<vector<pair<double, int> > > iqtreeRateParser(const string &filename){
        vector<vector<pair<double, int> > > result;
        string line;
        ifstream fin(filename);
        while (getline(fin, line)){
            if (line[0] == '#') continue;
            if (line.substr(0,4) == "Site"){
                result.emplace_back();
                double total = 0;
                while (getline(fin, line)){
                    if (line[0] == '#') break;
                    stringstream ss(line);
                    double r;
                    int i;
                    ss >> i >> r;
                    result.back().push_back({r, i - 1});
                    total += r;
                }
                for (auto& e: result.back()) e.first /= total / result.back().size();
            }
        }
        return result;
    }
}

struct TreeTokenizer{
    static const string KEYWORDS, SPACES;
    
    const string TEXT;
    int i = 0, j;

    TreeTokenizer(const string& TEXT): TEXT(TEXT){}

    string operator()(){
        string res = preview();
        i = j;
        return res;
    }

    string preview(){
        while (i < TEXT.size() && SPACES.find(TEXT[i]) != string::npos) i++;
        if (i == TEXT.size()) return "";
        string res;
        if (KEYWORDS.find(TEXT[i]) != string::npos){
            j = i+1;
            return TEXT.substr(i, 1);
        }
        else{
            bool singleQ = false, doubleQ = false;
            for (j = i; j < TEXT.size(); j++){
                if (TEXT[j] == '\'') singleQ = !singleQ;
                if (TEXT[j] == '\"') doubleQ = !doubleQ;
                if (!singleQ && !doubleQ && KEYWORDS.find(TEXT[j]) != string::npos) break;
            }
            int k = j - 1;
            while (k >= i && SPACES.find(TEXT[k]) != string::npos) k--;
            return TEXT.substr(i, k-i+1);
        }
    }
};
const string TreeTokenizer::KEYWORDS = "(),:;", TreeTokenizer::SPACES = " \t\n\r";

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
        for (int i = 0; i < s.size(); i++){
            y.s.push_back(s[i] | x[i]);
        }
        y.p = p + x.p;
        return y;
    }

    AlignmentHot operator-(const AlignmentHot &x) const{
        AlignmentHot y;
        for (int i = 0; i < s.size(); i++){
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
        for (int i = 1; i < s.size(); i++)
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
        int n = a.nTaxa(), L = a[0].nSites();
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                /*int x = 0, y = 0;
                for (int k = 0; k < L; k++){
                    if ((a[i][k] && b[j][k]) || (a[j][k] && b[i][k])) x++;
                    if ((a[i][k] || b[i][k]) && (a[j][k] || b[j][k])) y++;
                }*/
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

struct Tree{
    struct Node
    {
        bool isleaf = false;
        int parent = -1, taxon = -1;
        string label;
        double length = 0;
    };
    
    vector<Node> nodes;
    int root = -1;
    
    Tree(){}

    Tree(const string& TEXT){
        TreeTokenizer tk(TEXT);
        root = readSubtree(tk);
    }

    Tree(const string& TEXT, const unordered_map<string, int>& name2id): Tree(TEXT){
        for (Node &node: nodes){
            if (name2id.count(node.label)) node.taxon = name2id.at(node.label);
        }
    }

    Tree(DistanceMatrix dm){
        int n = dm.nTaxa(), m = n;
        if (n == 1){
            nodes.emplace_back();
            nodes[0].taxon = 0;
            root = 0;
            return;
        }
        vector<vector<double> > D = dm.d;
        vector<int> nodePos;
        for (int i = 0; i < n; i++){
            nodes.emplace_back();
            nodes.back().taxon = i;
            nodePos.push_back(i);
        }
        for (; n > 2; n--, m++){
            vector<double> S(n);
            double minq = 1e9, x = 0, y = 0;
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    S[i] += D[i][j];
            for (int i = 0; i < n; i++){
                for (int j = i + 1; j < n; j++){
                    double q = (n - 2) * D[i][j] - S[i] - S[j];
                    if (q < minq){ minq = q; x = i; y = j; }
                }
            }
            nodes.emplace_back();
            nodes[nodePos[x]].parent = m;
            nodes[nodePos[x]].length = max((D[x][y] + (S[x] - S[y]) / (n - 2)) / 2, 0.0);
            nodes[nodePos[y]].parent = m;
            nodes[nodePos[y]].length = max((D[x][y] + (S[y] - S[x]) / (n - 2)) / 2, 0.0);
            for (int i = 0; i < n; i++)
                if (i != x && i != y) D[i][x] = D[x][i] = (D[x][i] + D[y][i] - D[x][y]) / 2;
            nodePos[x] = m; nodePos[y] = -1; 
            vector<vector<double> > newD(n - 1);
            vector<int> newNodePos, translate;
            for (int i = 0; i < n; i++){
                if (i != y) {
                    newNodePos.push_back(nodePos[i]);
                    translate.push_back(i);
                }
            }
            for (int i = 0; i < n - 1; i++)
                for (int j = 0; j < n - 1; j++)
                    newD[i].push_back(D[translate[i]][translate[j]]);
            nodePos = newNodePos;
            D = newD;
        }
        nodes.emplace_back();
        nodes[nodePos[0]].parent = m;
        nodes[nodePos[0]].length = max(D[0][1], 0.0);
        nodes[nodePos[1]].parent = m;
        nodes[nodePos[1]].length = 0;
        root = m;
    }

    Tree(DistanceMatrix dm, const vector<int> &renumber): Tree(dm){
        for (Node &node: nodes){
            if (node.taxon != -1) node.taxon = renumber[node.taxon];
        }
    }

    int readSubtree(TreeTokenizer& tk){
        int v = nodes.size();
        nodes.emplace_back();
        string s = tk();
        if (s != "(") {
            nodes[v].isleaf = true;
            nodes[v].label = SeqUtils::PARSE_LEAFNAME(s);
            s = tk.preview();
            if (s == ":"){
                tk();
                s = tk();
                nodes[v].length = SeqUtils::from_string(s);
            }
        }
        else{
            do{
                int c = readSubtree(tk);
                nodes[c].parent = v;
            } while(tk() != ")");
            s = tk.preview();
            if (TreeTokenizer::KEYWORDS.find(s[0]) == string::npos) {
                nodes[v].label = tk();
                s = tk.preview();
            }
            if (s == ":"){
                tk();
                s = tk();
                nodes[v].length = SeqUtils::from_string(s);
            }
        }
        return v;
    }

    void diskCoveringNBinRecursion(int v, vector<int> &curbin, vector<vector<int> > &bins,
            const unordered_set<int> &breakpoints, const unordered_map<int, vector<int> > &children){
        if (breakpoints.count(v)){
            vector<int> newbin;
            if (nodes[v].isleaf){
                newbin.push_back(nodes[v].taxon);
            }
            else{
                for (int c: children.at(v)) diskCoveringNBinRecursion(c, newbin, bins, breakpoints, children);
            }
            bins.push_back(newbin);
        }
        else{
            if (nodes[v].isleaf){
                curbin.push_back(nodes[v].taxon);
            }
            else{
                for (int c: children.at(v)) diskCoveringNBinRecursion(c, curbin, bins, breakpoints, children);
            }
        }
    }

    vector<vector<int> > diskCoveringNBin(int nBin, mt19937 &eng){
        unordered_map<int, vector<int> > children;
        vector<int> branches;
        unordered_set<int> breakpoints;
        vector<vector<int> > bins;
        vector<int> lastBin;
        for (int i = 0; i < nodes.size(); i++){
            if (i != root) children[nodes[i].parent].push_back(i);
        }
        for (auto const &e: children){
            if (e.first == root && e.second.size() == 2){
                branches.push_back(e.second[0]);
            }
            else{
                for (int i: e.second) branches.push_back(i);
            }
        }
        shuffle(branches.begin(), branches.end(), eng);
        branches.resize(nBin - 1);
        for (int i: branches) breakpoints.insert(i);
        diskCoveringNBinRecursion(root, lastBin, bins, breakpoints, children);
        bins.push_back(lastBin);
        return bins;
    }

    int diskCoveringRecursion(int id, vector<vector<int> > &result, vector<int> &cluster,
            mt19937 &eng, double rate, uniform_real_distribution<double> &dist) const{
        if (cluster[id] != -1) return cluster[id];
        if (nodes[id].parent == -1 || dist(eng) > exp(-nodes[id].length * rate)) {
            cluster[id] = result.size();
            result.emplace_back();
        }
        else cluster[id] = diskCoveringRecursion(nodes[id].parent, result, cluster, eng, rate, dist);
        if (nodes[id].taxon != -1) result[cluster[id]].push_back(nodes[id].taxon);
        return cluster[id];
    }

    vector<vector<int> > diskCovering(mt19937 &eng, double rate = 1) const{
        vector<vector<int> > result;
        uniform_real_distribution<double> dist(0.0, 1.0);
        vector<int> cluster(nodes.size(), -1);
        for (int id = 0; id < nodes.size(); id++){
            if (nodes[id].taxon != -1) diskCoveringRecursion(id, result, cluster, eng, rate, dist);
        }
        return result;
    }

    vector<vector<int> > diskCovering(mt19937 &eng, double rate, int replicate) const{
        vector<vector<int> > result;
        for (int r = 0; r < replicate; r++){
            const vector<vector<int> > dc = diskCovering(eng, rate);
            for (const vector<int> &group: dc) {
                result.push_back(group);
            }
        }
        return result;
    }

    string subtree2string(int v, const vector<vector<int> > &children) const{
        string str = "(";
        bool first = true;
        if (children[v].size() == 0){
            return to_string(nodes[v].taxon) + ":" + to_string(nodes[v].length);
        }
        for (int c: children[v]){
            if (first) first = false;
            else str += ",";
            str += subtree2string(c, children);
        }
        return str + "):" + to_string(nodes[v].length);
    }

    operator string() const{
        vector<vector<int> > children(nodes.size());
        for (int i = 0; i < nodes.size(); i++){
            int p = nodes[i].parent;
            if (p != -1) children[p].push_back(i); 
        }
        return subtree2string(root, children) + ";";
    }

    string subtree2stringWithNames(int v, const vector<vector<int> > &children, const vector<string> &names, bool keepLength) const{
        string str = "(";
        bool first = true;
        if (children[v].size() == 0){
            if (!keepLength) return names[nodes[v].taxon];
            else return names[nodes[v].taxon] + ":" + to_string(nodes[v].length);
        }
        for (int c: children[v]){
            if (first) first = false;
            else str += ",";
            str += subtree2stringWithNames(c, children, names, keepLength);
        }
        if (!keepLength) return str + ")";
        else return str + "):" + to_string(nodes[v].length);
    }
    
    string stringWithNames(const vector<string> &names, bool keepLength = true) const{
        vector<vector<int> > children(nodes.size());
        for (int i = 0; i < nodes.size(); i++){
            int p = nodes[i].parent;
            if (p != -1) children[p].push_back(i); 
        }
        return subtree2stringWithNames(root, children, names, keepLength) + ";";
    }

    void minDistRecursion(int v, const vector<vector<int> > &children, vector<double> &minDist) const{
        if (children[v].size() != 2) {
            minDist[v] = nodes[v].length;
            return;
        }
        double temp = 1e6;
        for (int c: children[v]){
            minDistRecursion(c, children, minDist);
            temp = min(temp, minDist[c]);
        }
        minDist[v] = temp + nodes[v].length;
    }

    void shortQuartetLengthsRecursion(int v, const vector<vector<int> > &children, vector<double> &minDist,
            vector<double> &result, double upLength) const{
        if (children[v].size() != 2) return;
        int cl = children[v][0], cr = children[v][1];
        if (children[cl].size() == 2) {
            int clcl = children[cl][0], clcr = children[cl][1];
            result.push_back(minDist[clcl] + minDist[clcr] + nodes[cl].length + minDist[cr] + upLength);
            shortQuartetLengthsRecursion(cl, children, minDist, result, upLength + nodes[cl].length);
        }
        if (children[cr].size() == 2) {
            int crcl = children[cr][0], crcr = children[cr][1];
            result.push_back(minDist[crcl] + minDist[crcr] + nodes[cr].length + minDist[cl] + upLength);
            shortQuartetLengthsRecursion(cr, children, minDist, result, upLength + nodes[cr].length);
        }
    }

    vector<double> shortQuartetLengths() const{
        vector<double> result, minDist(nodes.size());
        vector<vector<int> > children(nodes.size());
        for (int i = 0; i < nodes.size(); i++){
            int p = nodes[i].parent;
            if (p != -1) children[p].push_back(i); 
        }
        minDistRecursion(root, children, minDist);
        if (children[root].size() == 2){
            int cl = children[root][0], cr = children[root][1];
            if (children[cl].size() == 2 && children[cr].size() == 2){
                int clcl = children[cl][0], clcr = children[cl][1];
                int crcl = children[cr][0], crcr = children[cr][1];
                result.push_back(minDist[clcl] + minDist[clcr] + minDist[crcl] + minDist[crcr]
                    + nodes[cl].length + nodes[cr].length);
            }
            shortQuartetLengthsRecursion(cl, children, minDist, result, nodes[cl].length + minDist[cr]);
            shortQuartetLengthsRecursion(cr, children, minDist, result, nodes[cr].length + minDist[cl]);
        }
        return result;
    }
};

vector<vector<int> > diskCovering(const AlignmentHot &a, const AlignmentHot &b,
        mt19937 &eng, double rate = 1, double threshold = 3, int replicates = 1){
    DistanceMatrix dm(a, b);
    vector<vector<int> > result;
    for (int r = 0; r < replicates; r++){
        vector<vector<int> > c = dm.cluster(eng, threshold);
        vector<DistanceMatrix> p = dm.partition(c);
        for (int i = 0; i < c.size(); i++){
            const Tree t(p[i]);
            //cerr << (string) t << endl;
            const vector<vector<int> > dc = t.diskCovering(eng, rate);
            for (const vector<int> &group: dc) {
                result.emplace_back();
                for (int id: group) result.back().push_back(c[i][id]);
            } 
        }
    }
    return result;
}

#endif