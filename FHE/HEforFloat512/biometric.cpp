#include "biometric.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <algorithm>
#include "dirent.h"

using namespace std;

double sciToDub(const string& str){ //https://www.oreilly.com/library/view/c-cookbook/0596007612/ch03s06.html
   stringstream ss(str);
   double d = 0;
   ss >> d;
   if (ss.fail()) {
      string s = "Unable to format ";
      s += str;
      s += " as a number!";
      throw (s);
   }
   return (d);
}

void readTemplate(int subjectID, string sample, double (&tempArr)[templateSize], int (&originalIDs)[nrSubj]){
    string paddedSubjectID = string(5 - to_string(originalIDs[subjectID]).length(), '0') + to_string(originalIDs[subjectID]);
    string s;
    if(subjectID < nrSubjFERET){
        s = "/Users/cecilie/Documents/Code/identification/HE_Features/feret_arcface512_float/" + paddedSubjectID + "_" + sample + ".af";
    } else s = "/Users/cecilie/Documents/Code/identification/HE_Features/frgc_arcface512_float/" + paddedSubjectID + sample + ".af";
    char * path = new char [s.size() + 1];
    strcpy(path, s.c_str());
    ifstream templateFS;
    templateFS.open(path);
    if (!templateFS) {
        cerr << "Unable to open file " + s;
        exit(1);
    }
    char buffer[30];
    string line;
    for(int i = 0; i < templateSize; i++){
        templateFS.getline(buffer, 30);
        line = buffer;
        tempArr[i] = sciToDub(line);
    }
    templateFS.close();
}

void encryptTemplate(Ciphertext<Element> &encTemplate, double (&plainTemplate)[templateSize], CryptoContext<Element> &cc, KeyPair<Element> &keyPair){
    vector<double> tmpVector (plainTemplate, plainTemplate + sizeof(plainTemplate) / sizeof(plainTemplate[0]));
    Plaintext tempP = cc->MakeCKKSPackedPlaintext(tmpVector);
    encTemplate = cc->Encrypt(keyPair.publicKey, tempP);
}

void encryptBatchedTemplate(Ciphertext<Element> &encTemplate, double (&plainTemplate)[8*templateSize], CryptoContext<Element> &cc, KeyPair<Element> &keyPair){
    vector<double> tmpVector (plainTemplate, plainTemplate + sizeof(plainTemplate) / sizeof(plainTemplate[0]));
    Plaintext tempP = cc->MakeCKKSPackedPlaintext(tmpVector);
    encTemplate = cc->Encrypt(keyPair.publicKey, tempP);
}

void setupEncDB(string (&subjects)[nrSubj], string (&samples)[nrSamp + 4], Ciphertext<Element> (&encDB)[nrSubj], int (&originalIDs)[nrSubj], CryptoContext<Element> &cc, KeyPair<Element> &keyPair){
    DIR *dir;
    struct dirent *ent;
    int i = 0;
    if ((dir = opendir ("/Users/cecilie/Documents/Code/identification/HE_Features/feret_arcface512_float/")) != NULL) {
        while ((ent = readdir (dir)) != NULL) {
            if (ent->d_name[0] != '.'){
                samples[i] = ent->d_name; //implicit char* to str
                i++;
            }
        }
        closedir (dir);
    }
    i = 0;
    if ((dir = opendir ("/Users/cecilie/Documents/Code/identification/HE_Features/frgc_arcface512_float/")) != NULL) {
        while ((ent = readdir (dir)) != NULL) {
            if (ent->d_name[0] != '.'){
                samples[nrSampFERET + i] = ent->d_name; //implicit char* to str
                i++;
            }
        }
        closedir (dir);
    }
    vector<string> v(samples, samples + sizeof(samples) / sizeof(samples[0]));
    sort(v.begin(), v.end());
    for(int j = 0; j < nrSamp; j++){
        samples[j] = v[j+4];
    }
    subjects[0] = samples[0];
    originalIDs[0] = stoi(samples[0].substr(0,5));
    int k = 1;
    for(int i = 1; i < nrSamp; i++){
        if(samples[i].substr(0,5) != samples[i-1].substr(0,5)){
            subjects[k] = samples[i];
            originalIDs[k] = stoi(samples[i].substr(0,5));
            k++;
        }
    }
    double tempTemplate[templateSize] = {};
    for(int i = 0; i < nrSubjFERET; i++){
        string s = "/Users/cecilie/Documents/Code/identification/HE_Features/feret_arcface512_float/" + subjects[i];
        char * path = new char [s.size() + 1];
        strcpy(path, s.c_str());
        ifstream templateFS;
        templateFS.open(path);
        if (!templateFS) {
            cerr << "Unable to open file " + s;
            exit(1);
        }
        char buffer[30];
        string line;
        for(int i = 0; i < templateSize; i++){
            templateFS.getline(buffer, 30);
            line = buffer;
            tempTemplate[i] = sciToDub(line);
        }
        encryptTemplate(encDB[i], tempTemplate, cc, keyPair);
        templateFS.close();
    }
    for(int i = 0; i < nrSubjFRGC; i++){
        string s = "/Users/cecilie/Documents/Code/identification/HE_Features/frgc_arcface512_float/" + subjects[nrSubjFERET + i];
        char * path = new char [s.size() + 1];
        strcpy(path, s.c_str());
        ifstream templateFS;
        templateFS.open(path);
        if (!templateFS) {
            cerr << "Unable to open file " + s;
            exit(1);
        }
        char buffer[30];
        string line;
        for(int i = 0; i < templateSize; i++){
            templateFS.getline(buffer, 30);
            line = buffer;
            tempTemplate[i] = sciToDub(line);
        }
        encryptTemplate(encDB[nrSubjFERET + i], tempTemplate, cc, keyPair);
        templateFS.close();
    }
}


void setupEncDBwBatching(string (&subjects)[nrSubj], string (&samples)[nrSamp + 4], Ciphertext<Element> (&encDB)[nrSubj], int (&originalIDs)[nrSubj], CryptoContext<Element> &cc, KeyPair<Element> &keyPair){
    DIR *dir;
    struct dirent *ent;
    int i = 0;
    if ((dir = opendir ("/Users/cecilie/Documents/Code/identification/HE_Features/feret_arcface512_float/")) != NULL) {
        while ((ent = readdir (dir)) != NULL) {
            if (ent->d_name[0] != '.'){
                samples[i] = ent->d_name; //implicit char* to str
                i++;
            }
        }
        closedir (dir);
    }
    i = 0;
    if ((dir = opendir ("/Users/cecilie/Documents/Code/identification/HE_Features/frgc_arcface512_float/")) != NULL) {
        while ((ent = readdir (dir)) != NULL) {
            if (ent->d_name[0] != '.'){
                samples[nrSampFERET + i] = ent->d_name; //implicit char* to str
                i++;
            }
        }
        closedir (dir);
    }
    vector<string> v(samples, samples + sizeof(samples) / sizeof(samples[0]));
    sort(v.begin(), v.end());
    for(int j = 0; j < nrSamp; j++){
        samples[j] = v[j+4];
    }
    subjects[0] = samples[0];
    originalIDs[0] = stoi(samples[0].substr(0,5));
    int k = 1;
    for(int i = 1; i < nrSamp; i++){
        if(samples[i].substr(0,5) != samples[i-1].substr(0,5)){
            subjects[k] = samples[i];
            originalIDs[k] = stoi(samples[i].substr(0,5));
            k++;
        }
    }
    double tempTemplate[templateSize] = {};
    for(int i = 0; i < 133; i++){
    //for(int i = 0; i < nrSubjFERET; i++){ //test for executions
        string s = "/Users/cecilie/Documents/Code/identification/HE_Features/feret_arcface512_float/" + subjects[i];
        char * path = new char [s.size() + 1];
        strcpy(path, s.c_str());
        ifstream templateFS;
        templateFS.open(path);
        if (!templateFS) {
            cerr << "Unable to open file " + s;
            exit(1);
        }
        char buffer[30];
        string line;
        for(int i = 0; i < templateSize; i++){
            templateFS.getline(buffer, 30);
            line = buffer;
            tempTemplate[i] = sciToDub(line);
        }
        encryptTemplate(encDB[i], tempTemplate, cc, keyPair);
        templateFS.close();
    }
    /* //uncomment to test for execution
    for(int i = 0; i < nrSubjFRGC; i++){
        string s = "/Users/cecilie/Documents/Code/identification/HE_Features/frgc_arcface512_float/" + subjects[nrSubjFERET + i];
        char * path = new char [s.size() + 1];
        strcpy(path, s.c_str());
        ifstream templateFS;
        templateFS.open(path);
        if (!templateFS) {
            cerr << "Unable to open file " + s;
            exit(1);
        }
        char buffer[30];
        string line;
        for(int i = 0; i < templateSize; i++){
            templateFS.getline(buffer, 30);
            line = buffer;
            tempTemplate[i] = sciToDub(line);
        }
        encryptTemplate(encDB[nrSubjFERET + i], tempTemplate, cc, keyPair);
        templateFS.close();
    }
     */
}