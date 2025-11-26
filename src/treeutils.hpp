#ifndef TREE_UTILS
#define TREE_UTILS

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

namespace TreeUtils{
    string PARSE_LEAFNAME(const string& TEXT){
        string s;
        for (int i = 0; i < TEXT.size(); i++){
            if (TEXT[i] != '\"' && TEXT[i] != '\'') s += TEXT[i];
        }
        return s;
    }

    double from_string(const string s){
    	return stof(s);
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

    Tree(TreeTokenizer &tk){
        root = readSubtree(tk);
        string s = tk();
        if (s != ";") {
            cerr << "Newick does not properly end with ';'.\n";
            exit(0);
        }
    }

    Tree(TreeTokenizer &tk, const unordered_map<string, int>& name2id): Tree(tk){
        for (Node &node: nodes){
            if (name2id.count(node.label)) node.taxon = name2id.at(node.label);
        }
    }

    Tree(const string& TEXT): Tree(){
        TreeTokenizer tk(TEXT);
        root = readSubtree(tk);
        string s = tk();
        if (s != ";") {
            cerr << "Newick does not properly end with ';'.\n";
            exit(0);
        }
    }

    Tree(const string& TEXT, const unordered_map<string, int>& name2id): Tree(TEXT){
        for (Node &node: nodes){
            if (name2id.count(node.label)) node.taxon = name2id.at(node.label);
        }
    }

    static vector<Tree> text2trees(const string& TEXT){
        TreeTokenizer tk(TEXT);
        vector<Tree> res;
        while (tk.preview() != ""){
            res.emplace_back(tk);
        }
        return res;
    }

    static vector<Tree> text2trees(const string& TEXT, const unordered_map<string, int>& name2id){
        TreeTokenizer tk(TEXT);
        vector<Tree> res;
        while (tk.preview() != ""){
            res.emplace_back(tk, name2id);
        }
        return res;
    }

    int readSubtree(TreeTokenizer& tk){
        int v = nodes.size();
        nodes.emplace_back();
        string s = tk();
        if (s != "(") {
            nodes[v].isleaf = true;
            nodes[v].label = TreeUtils::PARSE_LEAFNAME(s);
            s = tk.preview();
            if (s == ":"){
                tk();
                s = tk();
                nodes[v].length = TreeUtils::from_string(s);
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
                nodes[v].length = TreeUtils::from_string(s);
            }
        }
        return v;
    }

    operator string() const{
        return newick();
    }

    string subtree2string(int v, const vector<vector<int> > &children, bool keepLength, bool keepInternalLabel, bool useTaxonID) const{
        string str = "";
        bool first = true;
        if (children[v].size() == 0){
            str = ((useTaxonID) ? to_string(nodes[v].taxon) : nodes[v].label);
            if (keepLength) str += string(":") + to_string(nodes[v].length);
            return str;
        }
        str = "(";
        for (int c: children[v]){
            if (first) first = false;
            else str += ",";
            str += subtree2string(c, children, keepLength, keepInternalLabel, useTaxonID);
        }
        str += ")";
        if (keepInternalLabel) str += nodes[v].label;
        if (keepLength) str += string(":") + to_string(nodes[v].length);
        return str;
    }
    
    string newick(bool keepLength = true, bool keepInternalLabel = true, bool useTaxonID = false) const{
        vector<vector<int> > children(nodes.size());
        for (int i = 0; i < nodes.size(); i++){
            int p = nodes[i].parent;
            if (p != -1) children[p].push_back(i); 
        }
        return subtree2string(root, children, keepLength, keepInternalLabel, useTaxonID) + ";";
    }
    
    string subtree2string(int v, const vector<vector<int> > &children, const vector<string> &names, bool keepLength, bool keepInternalLabel) const{
        string str = "";
        bool first = true;
        if (children[v].size() == 0){
            if (!keepLength) return names[nodes[v].taxon];
            else return names[nodes[v].taxon] + ":" + to_string(nodes[v].length);
        }
        str = "(";
        for (int c: children[v]){
            if (first) first = false;
            else str += ",";
            str += subtree2string(c, children, names, keepLength, keepInternalLabel);
        }
        str += ")";
        if (keepInternalLabel) str += nodes[v].label;
        if (keepLength) str += string(":") + to_string(nodes[v].length);
        return str;
    }
    
    string newick(const vector<string> &names, bool keepLength = true, bool keepInternalLabel = true) const{
        vector<vector<int> > children(nodes.size());
        for (int i = 0; i < nodes.size(); i++){
            int p = nodes[i].parent;
            if (p != -1) children[p].push_back(i); 
        }
        return subtree2string(root, children, names, keepLength, keepInternalLabel) + ";";
    }
};

struct BinaryTree{
    struct Node{
        int parent = -1, taxon = -1, left = -1, right = -1;
        string label;
        double length = 0;
        bool ghost = false;
    };

    vector<Node> nodes;
    int root;
    
    int& parent(int i){
        return nodes[i].parent;
    }

    int& taxon(int i){
        return nodes[i].taxon;
    } 

    int& left(int i){
        return nodes[i].left;
    }

    int& right(int i){
        return nodes[i].right;
    }

    string& label(int i){
        return nodes[i].label;
    }

    double& length(int i){
        return nodes[i].length;
    }

    bool& ghost(int i){
        return nodes[i].ghost;
    }

    BinaryTree(){}

    BinaryTree(const Tree& tree, bool resolveNonbinary = true){
        vector<vector<int> > children(tree.nodes.size());
        for (int i = 0; i < tree.nodes.size(); i++){
            int p = tree.nodes[i].parent;
            if (p != -1) children[p].push_back(i); 
        }
        root = build(tree, children, tree.root, resolveNonbinary);
    }

    int build(const Tree& tree, const vector<vector<int> > &children, int tNode, bool resolveNonbinary){
        int bNode = nodes.size();
        nodes.emplace_back();
        if (!tree.nodes.at(tNode).isleaf){
            int bL = build(tree, children, children[tNode][0], resolveNonbinary);
            int bR = build(tree, children, children[tNode][1], resolveNonbinary);
            parent(bL) = bNode; left(bNode) = bL;
            parent(bR) = bNode; right(bNode) = bR;
            for (int i = 2; i < children[tNode].size(); i++){
                if (!resolveNonbinary) {
                    cerr << "Tree not binary!\n";
                    exit(0);
                }
                bL = bNode;
                ghost(bL) = true;
                bNode = nodes.size();
                nodes.emplace_back();
                bR = build(tree, children, children[tNode][i], resolveNonbinary);
                parent(bL) = bNode; left(bNode) = bL;
                parent(bR) = bNode; right(bNode) = bR;
            }
        }
        taxon(bNode) = tree.nodes.at(tNode).taxon;
        label(bNode) = tree.nodes.at(tNode).label;
        length(bNode) = tree.nodes.at(tNode).length;
        return bNode;
    }

    static vector<BinaryTree> text2trees(const string& TEXT){
        vector<Tree> trees = Tree::text2trees(TEXT);
        vector<BinaryTree> res;
        for (const Tree &t: trees) res.emplace_back(t);
        return res;
    }

    static vector<BinaryTree> text2trees(const string& TEXT, const unordered_map<string, int>& name2id){
        vector<Tree> trees = Tree::text2trees(TEXT, name2id);
        vector<BinaryTree> res;
        for (const Tree &t: trees) res.emplace_back(t);
        return res;
    }

    operator string() const{
        return newick(root);
    }

    string newick(bool keepLength = true, bool keepInternalLabel = true, bool useTaxonID = false) const{
        return subtree2string(root, keepLength, keepInternalLabel, useTaxonID) + ";";
    }

    string newick(const vector<string> &names, bool keepLength = true, bool keepInternalLabel = true) const{
        return subtree2string(root, names, keepLength, keepInternalLabel) + ";";
    }

private:
    string subtree2string(int v, bool keepLength, bool keepInternalLabel, bool useTaxonID) const{
        string str = "";
        bool first = true;
        if (nodes[v].left == -1){
            str = ((useTaxonID) ? to_string(nodes[v].taxon) : nodes[v].label);
            if (keepLength) str += string(":") + to_string(nodes[v].length);
            return str;
        }
        str = "(";
        str += subtree2string(nodes[v].left, keepLength, keepInternalLabel, useTaxonID);
        str += ",";
        str += subtree2string(nodes[v].right, keepLength, keepInternalLabel, useTaxonID);
        str += ")";
        if (keepInternalLabel) str += nodes[v].label;
        if (keepLength) str += string(":") + to_string(nodes[v].length);
        return str;
    }

    string subtree2string(int v, const vector<string> &names, bool keepLength, bool keepInternalLabel) const{
        string str = "";
        bool first = true;
        if (nodes[v].left == -1){
            str = names.at(nodes[v].taxon);
            if (keepLength) str += string(":") + to_string(nodes[v].length);
            return str;
        }
        str = "(";
        str += subtree2string(nodes[v].left, names, keepLength, keepInternalLabel);
        str += ",";
        str += subtree2string(nodes[v].right, names, keepLength, keepInternalLabel);
        str += ")";
        if (keepInternalLabel) str += nodes[v].label;
        if (keepLength) str += string(":") + to_string(nodes[v].length);
        return str;
    }
};
#endif