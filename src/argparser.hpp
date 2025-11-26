#ifndef ARG_PARSER
#define ARG_PARSER

#include<string>
#include<unordered_set>
#include<unordered_map>
#include<functional>
#include<iostream>

using namespace std;

#include "tutorial.hpp"

class ArgParser{
    unordered_set<string> longName;
    unordered_map<char,string> shortName;
    unordered_map<string,string> argType;
    unordered_map<string,string> argDescription;
    unordered_map<string,bool> argPriority;
    unordered_map<string,int> intValue;
    unordered_map<string,double> doubleValue;
    unordered_map<string,string> stringValue;
    unordered_map<string,function<void()> > func;
    string programName = "ASTER";
    string fullName = "Accurate Species Tree EstimatoR";

    void addArg(char sName, string lName, string description, bool priority = false){
        if (sName == 0) longName.insert(lName);
        else shortName[sName] = lName;
        argDescription[lName] = description;
        argPriority[lName] = priority;
    }

public:
    void addIntArg(char sName, string lName, int defaultValue, string description, bool priority = false){
        addArg(sName, lName, description, priority);
        argType[lName] = "Integer";
        intValue[lName] = defaultValue;
    }

    void addDoubleArg(char sName, string lName, double defaultValue, string description, bool priority = false){
        addArg(sName, lName, description, priority);
        argType[lName] = "Float";
        doubleValue[lName] = defaultValue;
    }

    void addStringArg(char sName, string lName, string defaultValue, string description, bool priority = false){
        addArg(sName, lName, description, priority);
        argType[lName] = "String";
        stringValue[lName] = defaultValue;
    }

    void addFlag(char sName, string lName, string description, function<void()> f, bool priority = false){
        addArg(sName, lName, description, priority);
        argType[lName] = "Preset";
        func[lName] = f;
    }

    int& getIntArg(string name){
        return intValue[name];
    }

    double& getDoubleArg(string name){
        return doubleValue[name];
    }

    string& getStringArg(string name){
        return stringValue[name];
    }

    void setProgramName(string name, string fname = "Accurate Species Tree EstimatoR"){
        programName = name;
        fullName = fname;
    }

    string getFullName(){
        return fullName;
    }

    void printHelp(string name){
        cerr << "\t--" << name << "\t" << argType[name] << "\t";
        if (argType[name].compare("Integer") == 0) cerr << intValue[name];
        if (argType[name].compare("Float") == 0) cerr << doubleValue[name];
        if (argType[name].compare("String") == 0) cerr << stringValue[name];
        cerr << "\t" << argDescription[name] << endl;
    }

    ArgParser(){
        addStringArg('i', "input", "<last argument>", "The input file");
        addFlag('h', "help", "Printing the help message", [&](){
            cerr << programName << " [ --ARG_NAME ARG_VALUE ... ] input_file" << endl;
            cerr << "Abbr\tArg\tType\tDefault\tDescription" << endl;
            for (auto e: shortName) if (argPriority[e.second] && argType[e.second].compare("Preset") == 0) {cerr << "-" << e.first; printHelp(e.second);}
            for (auto e: longName) if (argPriority[e] && argType[e].compare("Preset") == 0) printHelp(e);
            for (auto e: shortName) if (argPriority[e.second] && argType[e.second].compare("Preset") != 0) {cerr << "-" << e.first; printHelp(e.second);}
            for (auto e: longName) if (argPriority[e] && argType[e].compare("Preset") != 0) printHelp(e);
            cerr << endl << "Advanced options:" << endl;
            for (auto e: shortName) if (!argPriority[e.second] && argType[e.second].compare("Preset") == 0) {cerr << "-" << e.first; printHelp(e.second);}
            for (auto e: longName) if (!argPriority[e] && argType[e].compare("Preset") == 0) printHelp(e);
            for (auto e: shortName) if (!argPriority[e.second] && argType[e.second].compare("Preset") != 0) {cerr << "-" << e.first; printHelp(e.second);}
            for (auto e: longName) if (!argPriority[e] && argType[e].compare("Preset") != 0) printHelp(e);
            cerr << endl << "Examples: " << endl;
            cerr << programName << " -o output_file input_file" << endl;
            cerr << programName << " -i input_file -o output_file" << endl;
            exit(0);
        }, true);
    }

    void parse(int argc, char** argv){
        bool inputLast = true;
        if (argc == 1) func["help"]();
        for (int i = 1; i < argc; i++){
            if (strcmp(argv[i], "-H") == 0) {MDGenerator temp(programName); func["help"]();}
            if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) func["help"]();
            if (strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "--input") == 0) inputLast = false;
        }
        if (inputLast){
            stringValue["input"] = argv[argc - 1];
            argc--;
        }
        for (int i = 1; i < argc; i++){
            if (argv[i][0] != '-') {cerr << "Error: Input parsing error around `" << argv[i] << "`." << endl; func["help"]();}
            string name;
            if (argv[i][1] == '-') {
                name = &argv[i][2];
                if (argType.count(name) == 0) {cerr << "Error: Unknown argument `" << argv[i] << "`." << endl; func["help"]();}
            }
            else {
                char c = argv[i][1];
                if (shortName.count(c) == 0) {cerr << "Error: Unknown argument `" << argv[i] << "`." << endl; func["help"]();}
                name = shortName[c];
            }
            if (argType[name].compare("Preset") == 0) func[name]();
            if (argType[name].compare("Integer") == 0) intValue[name] = stoi(argv[++i]);
            if (argType[name].compare("Float") == 0) doubleValue[name] = stod(argv[++i]);
            if (argType[name].compare("String") == 0) stringValue[name] = argv[++i];
        }
        if (stringValue["input"].compare("<last argument>") == 0) {cerr << "Error: Missing input file." << endl; func["help"]();}
    }
} ARG;

#endif